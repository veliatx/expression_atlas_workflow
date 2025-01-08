#!/usr/bin/env nextflow

/*
* RNA-seq Download Pipeline
* Author: mitch@veliatx.com
* 
* This pipeline downloads RNA-seq data from either ENA (via FTP) or SRA (via AWS S3)
* and optionally subsamples the resulting FASTQ files.
*
* Key features:
* - Supports both single-end and paired-end reads
* - Downloads from ENA using wget or from SRA using AWS S3
* - Optional subsampling of FASTQ files using seqtk
* - Configurable error handling and resource allocation
*/

nextflow.enable.dsl = 2

/*
Download from ENA using FTP.
*/
process wget_ena {

    publishDir params.output_directory, mode: 'copy', overwrite: true, pattern: '*.gz'

    container params.container_fastqdownload

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

    publishDir params.output_directory, 
        mode: 'copy', 
        overwrite: true, 
        pattern: '*.gz', 
        enabled: !params.subsample

    container params.container_fastqdownload

    maxForks params.sra_to_fastq_forks

    input:
    tuple val(SRR_ID), val(SRX_ID), path(SRA), val(SINGLE_END)
    path(FASTA_REFERENCE)

    output:
    tuple val(SRR_ID), val(SRX_ID), val(SINGLE_END), path('*.fastq*')

    script:
    """
    #!/bin/bash

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
        for f in *.fastq
        do
            mv \\
                \${f} \\
                ${SRX_ID}_\${f}
            if [ '${params.subsample}' = 'false' ]
            then
                pigz -c -p ${params.pigz_threads} ${SRX_ID}_\${f} > ${SRX_ID}_\${f}.gz
            fi
        done
    else
        mv \\
            ${SRR_ID}.fastq \\
            ${SRX_ID}_${SRR_ID}.fastq
        if [ '${params.subsample}' = 'false' ]
        then
            pigz -c -p ${params.pigz_threads} ${SRX_ID}_${SRR_ID}.fastq > ${SRX_ID}_${SRR_ID}.fastq.gz
        fi
    fi

    if [ '${params.subsample}' = 'true' ]
    then
        if [ '${SINGLE_END}' = 'false' ]
        then
            for f in *.fastq
            do
                mv \\
                    \${f} \\
                    \${f%.*}.presample.fastq
            done
        else
            mv ${SRX_ID}_${SRR_ID}.fastq ${SRX_ID}_${SRR_ID}.presample.fastq
        fi
    fi
    """
}

/*
Subsample a fastq with seqtk to specified number reads. 
*/
process subsample_fastq {

    container params.container_fastqdownload

    publishDir params.output_directory, mode: 'copy', overwrite: true, pattern: '*.gz'

    maxForks params.subsample_fastq_forks

    input:
    tuple val(SRR_ID), val(SRX_ID), val(SINGLE_END), path(FASTQS), val(READ_COUNT)

    output:
    tuple val(SRR_ID), val(SRX_ID), path('*.fastq.gz')
    
    script:
    """
    #!/bin/sh

    if [ '${params.sampling_cutoff < READ_COUNT}' = 'true' ]
    then
        if [ '${SINGLE_END}' = 'true' ]
        then
            seqtk \\
                sample \\
                -s${params.sampling_seed} \\
                ${SRX_ID}_${SRR_ID}.presample.fastq \\
                ${params.sampling_depth} \\
                | \\
                pigz -c -p ${params.pigz_threads} \\
                    > ${SRX_ID}_${SRR_ID}.fastq.gz
        else
            for f in *.fastq
            do
                seqtk \\
                    sample \\
                    -s${params.sampling_seed} \\
                    \${f} \\
                    ${params.sampling_depth} \\
                    | \\
                    pigz -c -p ${params.pigz_threads} \\
                        > \${f%.*}.fastq.gz
            done
        fi
    else
        if [ '${SINGLE_END}' = 'true' ]
        then
            mv ${SRX_ID}_${SRR_ID}.presample.fastq ${SRX_ID}_${SRR_ID}.fastq
            pigz -c -p ${params.pigz_threads} \\
                ${SRX_ID}_${SRR_ID}.fastq \\
                > ${SRX_ID}_${SRR_ID}.fastq.gz
        else
           for f in *.fastq
            do
                mv \\
                    \${f} \\
                    \${f%%.*}.fastq
                pigz -c -p ${params.pigz_threads} \\
                    \${f%%.*}.fastq \\
                    > \${f%%.*}.fastq.gz
            done
        fi
    fi
    """
}

def help() {
    log.info """
    ===========================================
    RNA-seq Download Pipeline
    ===========================================
    
    Required Arguments:
        --samplesheet         Path to fetchngs metadata samplesheet
        --output_directory    Path to output directory
    
    Optional Arguments:
        --fasta_reference    Path to reference FASTA (for fastq-dump)
        --wget              Use wget to download from ENA instead of SRA
        --subsample         Enable FASTQ subsampling
        
    Subsampling Options:
        --sampling_seed     Random seed for reproducibility
        --sampling_cutoff   Minimum read count threshold for subsampling
        --sampling_depth    Target depth to subsample to
        
    Resource Control:
        --wget_ena_forks    Max parallel downloads from ENA
        --s3_sra_cp_forks   Max parallel downloads from SRA
        --pigz_threads      Number of threads for compression
    """
}

workflow {

    // Input validation
    if ( !params.samplesheet || !params.output_directory ) {
        help()
        exit 1
    }

    // Prepare reference FASTA if specified
    fa = params.fasta_reference ? file(params.fasta_reference) : file('nf_download/test/empty.txt')
    
    // Split workflow based on download method
    if (!params.wget) {
        // SRA download workflow
        Channel.fromPath(params.samplesheet)
            .splitCsv(header: true, sep: '\t')
            .map { row -> [
                row.run_accession,
                row.experiment_accession,
                row.single_end.toLowerCase()
            ]}
            .set { fastqs }

        // Main download processes
        s3_sra_cp(fastqs)
        sra_to_fastq(s3_sra_cp.out, fa)

        // Optional subsampling workflow
        if (params.subsample) {
            Channel.fromPath ( params.samplesheet )
                .splitCsv ( header: true, sep: '\t' )
                .map ( row -> [ row.run_accession, row.read_count as int ]) 
                .set { fastq_depths }
            sra_to_fastq.out
                .combine( fastq_depths, by: 0 )
                .set { downloaded_w_depths  }
            subsample_fastq( downloaded_w_depths )    
        }
    } else {
        // ENA download workflow
        Channel.fromPath ( params.samplesheet )
            .splitCsv ( header: true, sep: '\t' )
            .map ( row -> [row.run_accession, row.experiment_accession, row.fastq_1, row.fastq_2, row.single_end.toLowerCase()])
            .set { fastqs }

        wget_ena ( fastqs )
    }
}
