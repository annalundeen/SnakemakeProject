import argparse #handle command line arguments
import os #file handling
import gzip #unzip genome
import urllib.request #download genome (wget wasn't working)
from Bio import SeqIO #parse genbank file

parser = argparse.ArgumentParser() #crate argument parser
parser.add_argument("--index", required=True) #output index
parser.add_argument("--report", required=True) #output report
args = parser.parse_args()

# Create directories
os.makedirs(os.path.dirname(args.index), exist_ok=True) #create index directory
os.makedirs(os.path.dirname(args.report), exist_ok=True) #create report directory

genome_gb = "output/GCF_000845245.1.gbff"
#genome download link
url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/245/GCF_000845245.1_ViralProj14559/GCF_000845245.1_ViralProj14559_genomic.gbff.gz"

#download and unzip genome
if not os.path.exists(genome_gb):
    print(f"Downloading {genome_gb}...")
    gz_path = genome_gb + ".gz"
    urllib.request.urlretrieve(url, gz_path)
    with gzip.open(gz_path, 'rb') as f_in:
        with open(genome_gb, 'wb') as f_out:
            f_out.write(f_in.read())
    os.remove(gz_path)

#parse through genome
genome_fasta = "output/HCMV_genome.fa"
with open(genome_fasta, "w") as out:
    for record in SeqIO.parse(genome_gb, "genbank"):
        out.write(f">{record.id}\n{record.seq}\n")

#extract coding sequences
cds_fasta = "output/HCMV_CDS.fa"
cds_count = 0
with open(cds_fasta, "w") as out:
    for record in SeqIO.parse(genome_gb, "genbank"): #parse through genome
        for feature in record.features:
            if feature.type == "CDS" and "protein_id" in feature.qualifiers: #make sure record has cds and protein id- if it does, extract things and write them out
                seq = feature.extract(record.seq)
                pid = feature.qualifiers["protein_id"][0]
                out.write(f">{pid}\n{seq}\n")
                cds_count += 1 #count cds genes

#write out report
with open(args.report, "w") as f:
    f.write(f"The HCMV genome (GCF_000845245.1) has {cds_count} CDS.\n\n")

#build kallisto index
os.system(f"kallisto index -i {args.index} {cds_fasta}")


