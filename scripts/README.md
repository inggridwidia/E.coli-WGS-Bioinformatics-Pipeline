# Pipeline Automation Scripts

This folder contains the complete shell scripts used to automate the bioinformatics pipeline for *Escherichia coli* strain N31 Whole Genome Sequencing (WGS) analysis.

## Workflow Overview

1. **`mission01_data_preparation.sh`**: Quality control and raw data trimming.
2. **`mission02_genome_mapping.sh`**: Reference genome indexing and sequence alignment.
3. **`mission03_convert_sort.sh`**: SAM to BAM conversion, sorting, and indexing.
4. **`mission04_variant_calling.sh`**: Variant calling and annotating.
5. **`mission05_virulence_profiling.sh`**: Consensus generation, Prokka genome annotation, and BLASTp screening against VFDB.

---

## How to Run the Scripts

### Prerequisites
Ensure that all required bioinformatics tools (`SRA Toolkit`, `fastqc`. `fastp`, `minimap2`, `samtools`, `freebayes`, `bcftools`, `snpEff`, `prokka`, and `blastp`) are installed and properly configured in your Conda environment.


### Execution
All scripts use relative paths. Please always run them from the **project root directory**, not from inside the `scripts/` folder itself, to ensure all file paths resolve correctly.

```bash
# 1. Clone the repository and navigate to the root directory
cd E.coli-WGS-Bioinformatics-Pipeline

# 2. Activate your bioinformatics environment
conda activate your_env

# 3. Make sure the scripts have execution permissions
chmod +x scripts/*.sh

# 4. Run the desired script from the root directory. For example, to execute the genome annotation and virulence factor profiling
./scripts/mission05_virulence_profiling.sh
```

Note: For mission05_virulence_profiling.sh, the script will automatically check, download, and construct the BLAST database from the Virulence Factor Database (VFDB) if it is not already present in the data/vfdb/ directory.
