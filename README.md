# SnakemakeProject

## Dependencies
Must have latest version of Python installed: https://www.python.org/downloads/
Must have latest version of R installed: https://www.r-project.org/

conda install -c bioconda biopython snakemake kallisto bowtie2 spades blast-plus

## Installation
git clone https://github.com/annalundeen/SnakemakeProject.git
cd SnakemakeProject

## Running the Pipeline
Data in SnakemakeProject/data is test data, meaning a subset of data to run quickly. The final output should be PipelineReport.txt

To run pipeline:
snakemake --cores 4 

