#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

/*
Download from ENA using FTP.
*/
process wget_ena {

    publishDir params.output_directory, mode: 'copy', overwrite: true, pattern: '*.gz'

    container params.container_sratools

    maxForks params.wget_ena_forks

    input:
    tuple val(SRR_ID), val(SRX_ID), val(FTP_1), val(FTP_2), val(SINGLE_END)

    output:
    tuple val(SRR_ID), val(SRX_ID), path('*.gz'), val(SINGLE_END)

    script:
    """
    #!/bin/sh

    if [ '${SINGLE_END}' = 'false' ]
    then
        wget ${FTP_1}
        wget ${FTP_2}
    else
        wget ${FTP_1}
    fi
    """

}

/*
Copies SRR id out of public S3 bucket to local.
*/
process s3_sra_cp {

    errorStrategy params.errorStrategy

    container params.container_aws

    maxForks params.s3_sra_cp_forks

    input:
    tuple val(SRR_ID), val(SRX_ID), val(SINGLE_END)

    output:
    tuple val(SRR_ID), val(SRX_ID), path('*.sra'), val(SINGLE_END)

    script:
    """
    yum -yq install aws-cli 

    aws \\
        s3 \\
        cp \\
        s3://sra-pub-run-odp/sra/${SRR_ID}/${SRR_ID} \\
        ${SRR_ID}.sra \\
        --no-sign-request
    """
}

/*
Converts local SRA file to fastq file. Splits if paired-end.
*/
process sra_to_fastq {

    errorStrategy params.errorStrategy

    publishDir params.output_directory, mode: 'copy', overwrite: true, pattern: '*.fastq*'

    container params.container_sratools

    maxForks params.sra_to_fastq_forks

    input:
    tuple val(SRR_ID), val(SRX_ID), path(SRA), val(SINGLE_END)
    path(FASTA_REFERENCE)

    output:
    tuple val(SRR_ID), val(SRX_ID), path('*.fastq*')

    script:
    """
    #!/bin/sh

    if [ '${SINGLE_END}' = 'false' ]
    then
        fastq-dump \\
            ./${SRR_ID}.sra \\
            --split-files
    else
        fastq-dump \\
            ./${SRR_ID}.sra
    fi

    if [ '${SINGLE_END}' = 'false' ]
    then
        mv \\
            ${SRR_ID}_1.fastq \\
            ${SRX_ID}_${SRR_ID}_1.fastq
        mv \\
            ${SRR_ID}_2.fastq \\
            ${SRX_ID}_${SRR_ID}_2.fastq
        gzip *_1.fastq
       	gzip *_2.fastq
    else
        mv \\
            ${SRR_ID}.fastq \\
            ${SRX_ID}_${SRR_ID}.fastq
	gzip *.fastq
    fi
    """
}

def help() {
    log.info """
        Usage:
        nextflow run nf_download --samplesheet <path_to_samplesheet> --outdir <path_to_outdirectory>

        Arguments:
            --samplesheet fetchngs metadata samplesheet.
            --output_directory path to output directory. 
            --fasta_reference whether to force a reference fasta into fastq-dump. Provide path.
            --wget download from ena instead of sra.
        """
}

workflow {

    if ( !params.samplesheet || !params.output_directory ) {
        help()
        exit 1
    }

    if ( params.fasta_reference ) {
        fa = file( params.fasta_reference )
    } else {
        fa = file( 'download_sra/test/empty.txt' )
    }
    
    if ( !params.wget ) {
        Channel.fromPath ( params.samplesheet )
            .splitCsv ( header: true, sep: '\t' )
            .map ( row -> [row.run_accession, row.experiment_accession, row.single_end.toLowerCase()] )
            .set { fastqs }
    
        s3_sra_cp ( fastqs )
        sra_to_fastq ( s3_sra_cp.out, fa )
    } else {
        Channel.fromPath ( params.samplesheet )
            .splitCsv ( header: true, sep: '\t' )
            .map ( row -> [row.run_accession, row.experiment_accession, row.fastq_1, row.fastq_2, row.single_end.toLowerCase()])
            .set { fastqs }

        wget_ena ( fastqs)
    }
}
