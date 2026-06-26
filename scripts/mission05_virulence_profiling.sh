#!/bin/bash

set -e

echo "MISSION 5: Genome Annotation & Virulence Factor Profiling"

mkdir -p data/vfdb
mkdir -p mission05
mkdir -p mission05/annotation
mkdir -p mission05/virulence

# 1. Download and prepare VFDB database
if [ ! -f "data/vfdb/VFDB_setA_pro.pin" ]; then
    echo "Database VFDB not found. Downloading protein database (set A)..."
    wget -P data/vfdb/ http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
    
    echo "Extracting database file..."
    gunzip data/vfdb/VFDB_setA_pro.fas.gz

    echo "Creating a BLAST database from a VFDB file..." 
    makeblastdb -in data/vfdb/VFDB_setA_pro.fas -dbtype prot -out data/vfdb/VFDB_setA_pro
else
    echo "Database VFDB already exists. Skipping download."
fi

# 2. Consensus genome generation
echo "Compressing and indexing VCF file..."
bgzip -c mission04/variants_filtered.vcf > mission05/variants_filtered.vcf.gz
bcftools index mission05/variants_filtered.vcf.gz

echo "Generating consensus FASTA for Strain N31..."
bcftools consensus -f mission02/e_coli_RS218_reference.fna mission05/variants_filtered.vcf.gz > data/strain_N31_consensus.fasta

# 3. Whole genome annotation using Prokka
echo "Running Prokka genome annotation..."
prokka --outdir mission05/annotation \
       --prefix strain_N31 \
       --genus Escherichia \
       --species coli \
       --kingdom Bacteria \
       --force \
       data/strain_N31_consensus.fasta

# 4. Virulence factor screening against VFDB
echo "Running BLASTp screening against VFDB..."
blastp -query mission05/annotation/strain_N31.faa \
       -db data/vfdb/VFDB_setA_pro \
       -out mission05/virulence/virulence_blast_results.tsv \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -evalue 1e-5 \
       -num_threads 4

echo "Pipeline Completed Successfully!"
echo "Results saved to: mission05/virulence/virulence_blast_results.tsv"
