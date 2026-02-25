# Snakefile
#import packages
import os #handle paths and files

#define all input fastq files to avoid hard-coding later
SAMPLES = glob_wildcards("data/{sample}_1.fastq").sample
#define kallisto index location
INDEX = "output/kallisto/HCMV_index.idx"
#define final output report
REPORT = "PipelineReport.txt"

#rule to run everything needed to produce report
rule all:
    input:
        REPORT, #pipeline considered complete once we produce report file
        "output/sleuth/sleuth_results.tsv",
        expand("output/spades/{sample}/contigs.fasta", sample=SAMPLES), #assemblies for each sample
        expand("output/reports/{sample}_blast_report.txt", sample=SAMPLES) #blast reports per sample

rule compile_report:
    input:
        index_rep = "output/reports/index_report.txt",
        bowtie_reps = expand("output/reports/{sample}_bowtie_report.txt", sample = SAMPLES),
        blast_reps = expand("output/reports/{sample}_blast_report.txt", sample = SAMPLES),
        sleuth_rep = "output/sleuth/sleuth_results.tsv"
    output:
        REPORT
    shell:
        #printing out what we want for report: index summary, section headers, sleuth, bowtie, blast
        """
        cat {input.index_rep} > {output}
        echo -e '\n--- Sleuth Differential Expression ---' >> {output}
        cat {input.sleuth_rep} >> {output}
        echo -e '\n--- Bowtie2 Mapping ---' >> {output}
        cat {input.bowtie_reps} >> {output}
        echo -e '\n--- BLAST Results ---' >> {output}
        cat {input.blast_reps} >> {output}
        
        """

#rule to build kallisto index
rule build_index:
    output:
        index = INDEX, #produced by running build_index script
        genome = "output/HCMV_genome.fa",
        cds = "output/HCMV_CDS.fa",
        report = "output/reports/index_report.txt"
    shell: #running mython script
        "python3 scripts/build_index.py --index {output.index} --report {output.report}"

#rule to run kallisto quantification for each sample
rule quant_TPM:
    input:
        #get fastq paths from earlier sample dictionary using wildcard sample names
        fastq1 = "data/{sample}_1.fastq",
        fastq2 = "data/{sample}_2.fastq",
        index = INDEX #need kallisto index we built earlier
    output:
        "output/kallisto/{sample}/abundance.tsv" #kallisto output
    params:
        outdir = "output/kallisto/{sample}" #directory for kallisto results #directory for kallisto results
    shell: 
        "python3 scripts/quant_tpm.py --fastq1 {input.fastq1} --fastq2 {input.fastq2} --index {input.index} --outdir {params.outdir}" #running python command via shell

rule sleuth:
    input:
        expand("output/kallisto/{sample}/abundance.tsv", sample = SAMPLES)
    output:
        "output/sleuth/sleuth_results.tsv"
    shell:
        "Rscript scripts/run_sleuth.R {output}"

rule build_bowtie_index:
    input:
        "output/HCMV_genome.fa"
    output:
        "output/bowtie/HCMV_bowtie_index.1.bt2"
    params:
        prefix = "output/bowtie/HCMV_bowtie_index"
    shell:
        """
        mkdir -p output/bowtie
        bowtie2-build {input} {params.prefix}
        """

#map reads to HCMV genome and keep only mapped reads
rule bowtie_filter:
    input:
        #once again get fastq paths from sample directory using wildcard names
        fastq1 = "data/{sample}_1.fastq",
        fastq2 = "data/{sample}_2.fastq",
        genome = "output/HCMV_genome.fa",
        index = "output/bowtie/HCMV_bowtie_index.1.bt2"
    output:
        #output mapped files to appropriate sample directory
        mapped1 = "output/bowtie/{sample}_mapped_1.fastq",
        mapped2 = "output/bowtie/{sample}_mapped_2.fastq",
        report = "output/reports/{sample}_bowtie_report.txt"
    params:
        #need bowtie index, sam, and bam
        index_prefix = "output/bowtie/HCMV_bowtie_index",
        sam = "output/bowtie/{sample}.sam",
        bam = "output/bowtie/{sample}.bam"
    shell: #running bowtie filtering script
        "python3 scripts/bowtie_filter.py --fastq1 {input.fastq1} --fastq2 {input.fastq2} --index {params.index_prefix} --sam {params.sam} --bam {params.bam} --out1 {output.mapped1} --out2 {output.mapped2} --report {output.report}"

#rule to build transcriptome assemblies from our filtered reads
rule assemble_spades:
    input:
        #grab mapped samples instead of raw samples
        fastq1 = "output/bowtie/{sample}_mapped_1.fastq",
        fastq2 = "output/bowtie/{sample}_mapped_2.fastq"
    output:
        #output contigs
        contigs = "output/spades/{sample}/contigs.fasta"
    params:
        outdir = "output/spades/{sample}" #specifying output directory
    shell: #shell command to run spades python file
        """
        mkdir -p {params.outdir}

        python3 scripts/assemble_spades.py --fastq1 {input.fastq1} --fastq2 {input.fastq2} --outdir {params.outdir}

        ls -lh {params.outdir}

        if [ ! -f "{output.contigs}" ]; then
            echo "ERROR: contigs.fasta missing"
            exit 1
        fi
        """

rule build_blast_db:
    output:
        db_file = "blast_db/betaherpesvirinae_nr.nhr",
        fasta = "blast_db/betaherpesvirinae.fasta"
    params:
        db_name = "blast_db/betaherpesvirinae_nr"
    shell:
        """
        mkdir -p blast_db
        datasets download virus genome taxon 10359 --include genome --filename blast_db/virus.zip
        unzip -o blast_db/virus.zip -d blast_db/
        find blast_db/ncbi_dataset/data -name "*.fna" -exec cat {{}} + > {output.fasta}
        makeblastdb -in {output.fasta} -dbtype nucl -out {params.db_name}
        """

#rule to blast the longest contig and write out the top hits
rule blast_longest:
    input:
        #specify contigs from spades assembly
        contigs = "output/spades/{sample}/contigs.fasta",
        db_check = "blast_db/betaherpesvirinae_nr.nhr"
    output:
        #create output file for each sample when its done
        done_flag = "output/spades/{sample}/_blast_done.txt",
        report = "output/reports/{sample}_blast_report.txt"
    params:
        blast_db = "blast_db/betaherpesvirinae_nr" #must exist before running
    shell: #running blast python script
        "python3 scripts/blast_longest.py --contigs {input.contigs} --sample {wildcards.sample} --blast_db {params.blast_db} --report {output.report} --done {output.done_flag}"

rule clean:
    shell:
        "rm -rf output/ PipelineReport.txt .snakemake/"

#run as snakemake --cores 1
#snakemake -p --cores 4 --latency-wait 60
#to run cleanup rule: snakemake -c 1 clean
#dryrun snakemake -n
