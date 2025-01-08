# Nextflow RNA-seq Download Pipeline

A quick pipeline for downloading RNA-seq data from SRA (via AWS S3) or ENA (via FTP) with optional subsampling capabilities.

## Features

- Supports both single-end and paired-end reads
- Flexible download sources:
  - SRA using AWS S3 (default, no credentials needed)
  - ENA using direct FTP downloads
- Optional FASTQ subsampling using seqtk
- Configurable parallel processing
- Automatic handling of corrupted SRA files
- Compression using pigz for improved performance

## Usage

### Basic Usage with SRA (Default)
`nextflow run nf_download --samplesheet <samplesheet_path> --output_directory <output_dir> -w <work_dir> -profile docker`

### Using ENA Instead of SRA
`nextflow run nf_download --samplesheet <samplesheet_path> --output_directory <output_dir> -w <work_dir> --wget -profile docker`

### AWS Batch Execution
`nextflow run nf_download --samplesheet <samplesheet_path> --output_directory <output_dir> -w <work_dir> -profile aws`

## Samplesheet Format

Tab-separated file with the following columns:
- `run_accession`: SRR ID
- `experiment_accession`: SRX ID
- `single_end`: true/false
- `fastq_1`: FTP URL (required for ENA)
- `fastq_2`: FTP URL (required for paired-end ENA)
- `read_count`: Number of reads (required for subsampling)

## Notes

- Pipeline uses 'ignore' errorStrategy for handling corrupted SRA files
- Low default fork values for S3/ENA downloads to avoid rate limiting
- AWS Batch users should increase `sra_to_fastq_forks` for better performance
- Uses pigz for parallel compression of output files 
