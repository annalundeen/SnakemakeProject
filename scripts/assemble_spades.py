# 5. Using the Bowtie2 output reads, use SPAdes to generate four assemblies, one for each sample. You should use a k-mer size of 127. 

#import packages
import argparse #take in arguments and work with snakemake
import os #run things on command line

#adapted from homework argparse script, defining comamand line arguments
parser = argparse.ArgumentParser() #define parser
parser.add_argument("--fastq1", required = True) #add fastq1 argument
parser.add_argument("--fastq2", required = True) #add fastq2 argument
parser.add_argument("--outdir", required = True) #add outdir argument
args = parser.parse_args()

#make sure we have corect ourdirectory
os.makedirs(args.outdir, exist_ok = True)

#skip empty files
size1 = os.path.getsize(args.fastq1)
size2 = os.path.getsize(args.fastq2)

if size1 == 0 or size2 == 0:
    open(f"{args.outdir}/contigs.fasta", "w").close()
    exit()

#run spades to generate assemblies with k=127
exit_code=os.system(
f"spades.py -1 {args.fastq1} -2 {args.fastq2} -o {args.outdir} -k 127"
)

# ignore warnings exit code
if not os.path.exists(f"{args.outdir}/contigs.fasta"):
    exit(1)
exit(0)

