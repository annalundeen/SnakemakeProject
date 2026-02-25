#import packages
import argparse #take in arguments and work with snakemake
import os #running things on command line

#adapted from homework argparse script, defining comamand line arguments
parser = argparse.ArgumentParser() #define parser
parser.add_argument("--fastq1", required = True) #add fastq1 argument
parser.add_argument("--fastq2", required = True) #add fastq2 argument
parser.add_argument("--index", required = True) #add index argument
parser.add_argument("--outdir", required = True) #add outdir argument
args = parser.parse_args()

#make sure we have correct output directory
os.makedirs(args.outdir, exist_ok=True)

#run kallisto quant tpm command via command line
os.system(f"kallisto quant -i {args.index} -o {args.outdir} -b 30 -t 4 {args.fastq1} {args.fastq2}")


