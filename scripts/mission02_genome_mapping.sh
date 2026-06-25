#!/bin/bash

set -e

echo "MISSION 2: Genome Mapping"

mkdir -p mission02

# 1. Download Reference Genome
if [ ! -f "mission02/e_coli_RS218_reference.fna" ]; then
    echo "Downloading Escherichia coli RS218 reference genome..."
    wget -P mission02/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/817/345/GCF_000817345.1_ASM81734v1/GCF_000817345.1_ASM81734v1_genomic.fna.gz
    
    echo "Extracting reference genome..."
    gunzip mission02/GCF_000817345.1_ASM81734v1_genomic.fna.gz
    
    echo "Renaming reference genome file..."
    mv mission02/GCF_000817345.1_ASM81734v1_genomic.fna mission02/e_coli_RS218_reference.fna
else
    echo "Reference genome file already exists. Skipping download."
fi

# 2. Indexing Reference Genome
echo "Indexing reference genome..."
if [ ! -f "mission02/e_coli_RS218.mmi" ]; then
    minimap2 -d mission02/e_coli_RS218.mmi mission02/e_coli_RS218_reference.fna
else
    echo "Reference index already exists. Skipping indexing."
fi

# 3. Mapping Reads using minimap2
echo "Mapping reads to reference genome..."
minimap2 -ax map-hifi mission02/e_coli_RS218.mmi mission01/SRR36640870_clean.fastq > mission02/alignment.sam

echo "Mission 2 Completed Successfully!"
