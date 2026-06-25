#!/bin/bash

set -e

echo "MISSION 1: Data Preparation"

mkdir -p mission01

# 1. Download dataset from NCBI SRA
if [ ! -f "mission01/SRR36640870.fastq" ]; then
    echo "Downloading data SRR36640870..."
    prefetch SRR36640870
    
    echo "Converting SRA to FASTQ format..." 
    fastq-dump --outdir mission01 SRR36640870/SRR36640870.sra
else
    echo "Data file already exists. Skipping download."
fi

# 2. Raw read quality assessment
echo "Running Initial FastQC..."
fastqc mission01/SRR36640870.fastq -o mission01/

echo "FastQC analysis complete. Reports are in: mission01/"

# 3. Trimming with fastp
echo "Running Quality Trimming with fastp..."
fastp -i mission01/SRR36640870.fastq \
      -o mission01/SRR36640870_clean.fastq \
      --qualified_quality_phred 20 \
      --length_required 1000 \
      --html mission01/fastp_report.html \
      --json mission01/fastp_report.json

echo "Mission 1 Completed Successfully!"
