# Downloads fastq from sra/ena
### Example:
`nextflow run nf_download --samplesheet <samplesheet_path> --output_directory <output_dir> -w <work_dir> -profile docker`
### Example aws:
`nextflow run nf_download --samplesheet <samplesheet_path> --output_directory <output_dir> -w <work_dir> -profile aws`
### Notes:
* errorStrategy set to 'ignore' for entire pipeline -> corrupted fastqs coming out of sra won't crash pipeline
* s3_to_fastq_forks sets maxForks directive in s3_sra_cp process, defaults to 1 so that we fly under the radar
* wget_ena_forks sets maxForks directive in wget_ena process, defaults to 1 for same reason as above
* sra_to_fastq_forks sets maxForks directive in the fastq-dump process, defaults to 5 but should be increased if using batch 
