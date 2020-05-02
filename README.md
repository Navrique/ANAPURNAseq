# Nextflow pipeline for RNAseq read filtering and read count

> #Bioinformatics #pipeline 

## Tasks

- [ ] Test the pipeline on the a dataset  in [Data/Sample_dataset](Data/Sample_dataset)
  - [ ] adapt the code to reduce space usage at the rRNA step
  - [ ] test on [Data/Sample_dataset](Data/Sample_dataset)
- [ ] [README](Bioinformatics/ReadCountPipeline/README.md)
  - [ ] Usage 
  - [ ] Step description
  - [ ] introduction
- [ ] publish the pipeline on Github

## Done

- [x] Create a pipeline using [NextFlow](Notes/Nextflow.md) The pipeline can be found into [/Bioinformatics/ReadCountPipeline](/Bioinformatics/ReadCountPipeline)
  - [x] download all files
  - [x] generate fastQC reports
  - [x] generate MultiQC report
  - [x] sort reads using Clumpify
  - [x] Remove adapters 
  - [x] Quality filtring
  - [x] Artifact filtering
  - [x] Contamination removal
  - [x] filtering for low complexity
  - [x] Remove rRNA
  - [x] QC Analysis