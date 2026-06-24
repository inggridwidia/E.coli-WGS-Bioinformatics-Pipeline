#!/bin/bash

echo "MISSION 1: Data Preparation"

# 1. Download dataset from NCBI SRA
if [ ! -f "SRR36640870.sra" ]; then
    echo "Downloading data SRR36640870..."
    prefetch SRR36640870
    
    echo "Converting SRA to FASTQ format..." 
    fastq-dump SRR36640870/SRR36640870.sra
else
    echo "Data file already exists. Skipping download."
fi

# 2. Raw read quality assessment
echo "Running Initial FastQC..."
fastqc SRR36640870.fastq

echo "FastQC analysis complete. Individual reports are in: mission1/"

# 3. Trimming with fastp
echo "Running Quality Trimming with fastp..."
fastp -i SRR36640870.fastq \
      -o SRR36640870_clean.fastq \
      --qualified_quality_phred 20 \
      --length_required 1000 \
      --html fastp_report.html \
      --json fastp_report.json

echo "Mission 1 Completed Successfully!"
