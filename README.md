# Expression Atlas Download/Processing Workflow

This repository contains workflows for processing and analyzing RNA-seq data for the Expression Atlas.

## Download/Secondary Analysis Workflow

Located in `processing_workflow/`, this pipeline handles raw data acquisition and initial processing:
* Utilizes Nextflow and Jupyter notebooks for orchestration
* Fetches experiment metadata using nf-core/fetch_ngs
* Downloads raw sequencing data from SRA/ENA using custom Nextflow workflow
* Processes RNA-seq data using nf-core/rnaseq pipeline, including:
  - Quality control (FastQC)
  - Read trimming
  - Alignment
  - Quantification
* The processing workflow handles intermediate data processing steps:
    - Quality assessment of aligned reads
    - Normalization of expression values
    - Sample-level QC metrics
    - Data formatting and standardization
    - Generation of processing reports

## Differential Expression (DE) Workflow

Located in `de_workflow/`, the DE workflow performs downstream analysis:

* Experimental metadata curation and validation
* Quality control visualization and filtering
* Differential expression analysis including:
  - Fold change calculations
  - Statistical significance testing
  - Multiple testing correction
* Generation of results matrices and visualization
* Export of analysis artifacts in standardized formats
  - Standardized format required for ingestion into expression_atlas_db

## Dependencies

* Nextflow
* nf-core/fetch_ngs
* nf-core/rnaseq
* python environment created from environment.yml

## Usage
1. Copy `processing_workflow/` and `de_workflow/` to a project folder named after an SRA, ENA, or internal project id. 
  * Internal projects require creation of manifests that mimic those created by nf-core/fetchngs. See `processing_workflow/nf_download/test` for examples. 
2. Configure and run through 00_run_align_pipeline.ipynb in `processing_workflow/`. 
  * Build the nf_download container, if not build previously:
  ```
  docker build -f processing_workflow/nf_download/docker/Dockerfile -t fastq_download
  ```
  * Run individual nextflow commands as part of notebook, or copy nextflow commands to screens and run independently.
3. Configure and run through notebooks in order in `de_workflow/`.
4. Sync entire experiment folder to location in S3 for ingestion by expression_atlas_db. 

## Notes:
* As currently configured, this workflow designed to be ran on large local workstation (>500gb ram, >60 cores), but could be easily adapted to run on AWS batch. 