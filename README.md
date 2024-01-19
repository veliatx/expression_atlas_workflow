# Expression atlas download/processing workflow
### Download/Secondary Analysis Workflow 
* ipynb and nextflow workflow in download_workflow.
* relies on nf-core fetch_ngs for fetching metadata, a custom nextflow workflow for downloading from sra/ena, and nf-core rna-seq for the bulk of the processing. 
### Tertiary Analysis Workflow 
* ipynbs for defining metadata, qc, running de.
