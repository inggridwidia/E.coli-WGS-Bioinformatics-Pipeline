#!/bin/bash

set -e

echo "MISSION 3: Convert & Sort the Data"

mkdir -p mission03

# 1. Converting SAM to BAM
echo "Converting SAM to compressed BAM format..."
samtools view -Sb mission02/alignment.sam > mission03/alignment.bam

# 2. Sorting BAM file
echo "Sorting BAM file by genomic coordinates..."
samtools sort mission03/alignment.bam -o mission03/alignment.sorted.bam

# 3. Indexing sorted BAM file
echo "Generating BAM index (.bai)..."
samtools index mission03/alignment.sorted.bam

# 4. Generate simple statistics
echo "Generating alignment statistics..."
samtools flagstat mission03/alignment.sorted.bam
samtools depth mission03/alignment.sorted.bam | head -n 10
samtools coverage mission03/alignment.sorted.bam | head -n 10

echo "Mission 3 Completed Successfully!"
