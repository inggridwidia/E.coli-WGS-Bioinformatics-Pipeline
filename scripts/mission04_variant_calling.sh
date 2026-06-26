#!/bin/bash

set -e

echo "MISSION 4: Variant Calling & Annotation"

mkdir -p mission04 

# 1. Variant calling
echo "Running Variant Calling via FreeBayes..."
freebayes -f mission02/e_coli_RS218_reference.fna \
          --min-mapping-quality 20 \
          --min-base-quality 20 \
          --min-coverage 10 \
          mission03/alignment.sorted.bam > mission04/variants_raw.vcf

echo "Filtering VCF using bcftools..."
bcftools filter -e 'QUAL < 30 || INFO/DP < 10 || AF < 0.1' \
          mission04/variants_raw.vcf \
          -o mission04/variants_filtered.vcf

echo "Normalizing VCF using bcftools..."
bcftools norm -f mission02/e_coli_RS218_reference.fna \
          -m - \
          mission04/variants_filtered.vcf \
          -o mission04/variants_norm.vcf

# 2. Annotating variants
echo "Setting up local snpEff configuration & database..."
echo "e_coli_RS218.genome : Escherichia coli RS218" > snpEff.config
mkdir -p data/e_coli_RS218/
cp mission02/e_coli_RS218_reference.fna data/e_coli_RS218/sequences.fa

if [ ! -f data/e_coli_RS218/genes.gff ]; then
    echo "Downloading official RS218 GFF file from NCBI FTP..."
    wget -O mission04/e_coli_RS218_reference.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/817/345/GCF_000817345.1_ASM81734v1/GCF_000817345.1_ASM81734v1_genomic.gff.gz
    gunzip -f mission04/e_coli_RS218_reference.gff.gz
    cp mission04/e_coli_RS218_reference.gff data/e_coli_RS218/genes.gff
fi

echo "Building the custom snpEff database..."
snpEff build -gff3 -c snpEff.config -v e_coli_RS218

echo "Annotating Variants using snpEff..."
snpEff ann -c snpEff.config \
           -v \
           -stats mission04/snpEff_summary.html \
           e_coli_RS218 \
           mission04/variants_norm.vcf > mission04/variants_annotated.vcf

echo "Pipeline Completed Successfully! Preview of the annotated VCF:"
grep "^#CHROM" mission04/variants_annotated.vcf && grep -v "^#" mission04/variants_annotated.vcf | head -n 5
