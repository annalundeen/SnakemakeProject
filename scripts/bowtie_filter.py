#import packages
import argparse #take in arguments and work with snakemake
import os #running things on command line

#adapted from homework argparse script, defining command line arguments
parser = argparse.ArgumentParser() #define parser
parser.add_argument("--fastq1", required = True) #add fastq1 argument
parser.add_argument("--fastq2", required = True) #add fastq2 argument
parser.add_argument("--index", required = True) #add index argument
parser.add_argument("--sam", required = True) #add sam argument
parser.add_argument("--bam", required = True) #add bam argument
parser.add_argument("--out1", required = True) #add out1 argument
parser.add_argument("--out2", required = True) #add out2 argument
parser.add_argument("--report", required = True) #add report argument
args = parser.parse_args()

os.makedirs(os.path.dirname(args.report), exist_ok = True)

#count reads before mapping
before = 0 #intiialize before as 0
with os.popen(f"cat {args.fastq1} | wc -l") as p:
    before = int(p.read().strip()) // 4

#map reads
os.system(f'bowtie2 --very-sensitive -x {args.index} -1 {args.fastq1} -2 {args.fastq2} -S {args.sam}')
#convert sam to bam
os.system(f'samtools view -bS {args.sam} > {args.bam} 2>/dev/null')
#output mapped reads as fastq file
os.system(f'samtools fastq -f 2 -1 {args.out1} -2 {args.out2} {args.bam} 2>/dev/null')

#count reads after mapping 
after = int(os.popen(f'cat {args.out1} | wc -l').read()) // 4
sample = os.path.basename(args.out1).split("_mapped")[0]

#write out results to report
with open(args.report, "a") as f:
    f.write(f"Sample {sample} had {before} read pairs before and {after} read pairs after Bowtie2 filtering.\n")

