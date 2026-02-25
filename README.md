# SnakemakeProject

## Dependencies
Must have latest version of Python installed: https://www.python.org/downloads/
Must have latest version of R installed: https://www.r-project.org/

conda install -c bioconda biopython snakemake kallisto bowtie2 spades blast-plus

## Installation
git clone https://github.com/annalundeen/SnakemakeProject.git

cd SnakemakeProject

## Data Download
Test data (only first 10,000 reads) is avaiable in SnakemakeProject/data. Fastq files were downloaded from NCBI using fasterq-dump as follows.

prefetch SRR5660030
prefetch SRR5660033
prefetch SRR5660044
prefetch SRR5660045

fasterq-dump --split-files SRR5660030 -O /SnakemakeProject/data
fasterq-dump --split-files SRR5660033 -O /SnakemakeProject/data
fasterq-dump --split-files SRR5660044 -O /SnakemakeProject/data
fasterq-dump --split-files SRR5660045 -O /SnakemakeProject/data


## Running the Pipeline
Data in SnakemakeProject/data is test data, meaning a subset of data to run quickly. The final output should be PipelineReport.txt

To run pipeline:

cd SnakemakeProject

snakemake --cores 4 

