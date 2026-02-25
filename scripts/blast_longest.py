#importing packages we need
import argparse #take in arguments and work with snakemake
import os #running things on command line
from Bio import SeqIO #parse through fasta files

#adapted from homework argparse script, defining comamand line arguments
parser = argparse.ArgumentParser() #define parser
parser.add_argument("--contigs", required = True) #add contigs argument
parser.add_argument("--sample", required = True) #add sample argument
parser.add_argument("--blast_db", required = True) #add blast database argument
parser.add_argument("--report", required = True) #add report argument
parser.add_argument("--done", required = True)
args = parser.parse_args()

#make sure report_dir exists, make it if we don't have it
report_dir = os.path.dirname(args.report)
if report_dir:
    os.makedirs(report_dir, exist_ok=True)
os.makedirs("output/blast", exist_ok = True)

#skip empty contigs
if os.path.getsize(args.contigs) == 0:
    with open(args.report, "w") as f:
        f.write(f"{args.sample}: No contigs available for BLAST\n")
    open(args.done, "w").close()
    exit()

#find longest contig
longest_seq = None
max_len = 0

#parse through fasta file and find record with longest contig (max length)
for record in SeqIO.parse(args.contigs, "fasta"):
    if len(record.seq) > max_len:
        max_len = len(record.seq)
        longest_seq = record

#check if any contigs were found
if longest_seq is None:
    with open(args.report, "w") as f:
        f.write(f"{args.sample}: No contigs found\n")
    open(args.done, "w").close()
    exit()

#save longest contig
query_fasta = f"output/blast/{args.sample}_longest.fasta"

with open(query_fasta, "w") as f:
    f.write(f">{longest_seq.id}\n{longest_seq.seq}\n")

#define blast output location
blast_out = f"output/blast/{args.sample}_blast.txt"

#run blast on command line, output data we want
cmd = f"""
blastn \
-query {query_fasta} \
-db {args.blast_db} \
-out {blast_out} \
-outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" \
-max_target_seqs 5
"""
#storing exist code here if blast fails
ret = os.system(cmd)

if ret != 0: #if blast returned non-zero exit code (error)
    raise RuntimeError("BLAST failed")

#write top 5 blast hits to pipelinereport.txt
with open(args.report, "w") as f_out:
    f_out.write(f"{args.sample} BLAST Hits:\n")
    f_out.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")

    if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
        with open(blast_out) as f_in:
            for i, line in enumerate(f_in):
                if i >= 5:
                    break
                f_out.write(line)
    else:
        f_out.write("No significant hits found.\n")

open(args.done, "w").close()
