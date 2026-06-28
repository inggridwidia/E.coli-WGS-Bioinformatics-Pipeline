# ==============================================================================
# MISSION 05: VIRULENCE FACTOR VISUALIZATION & ANALYSIS
# Sample: E. coli Strain N31 (PacBio Long-Reads)
# Data Source: BLASTp against VFDB Set A
# ==============================================================================

# Load the required libraries
library(tidyverse)
library(ggplot2)

# PART 1: Import BLASTp Data File
# Define the relative path to the BLASTp output file
file_path <- "results/virulence/virulence_blast_results.tsv"

# Read the file and assign appropriate column names
blast_data <- read_tsv(
  file = file_path,
  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", 
                "gapopen", "qstart", "qend", "sstart", "send", 
                "evalue", "bitscore")
)

# Inspect the structure of the imported dataset
glimpse(blast_data)


# PART 2: Filter Valid Pathogenic Features
# Apply strict filtering thresholds
filtered_virulence <- blast_data %>% 
  filter(pident >= 80 & evalue <= 1e-5)

# Check the number of rows that passed the strict filter
nrow(filtered_virulence)

# Extract and standardise VFDB identifiers for downstream analysis
cleaned_virulence <- filtered_virulence %>% 
  mutate(vfid = str_extract(sseqid, "^VFG[0-9]+"))

# Preview the top 5 rows of the cleaned dataset
head(cleaned_virulence, 5)


# PART 3: Calculate Dominant Virulence Factors
# Compute the frequency of each VFDB identifier
virulence_counts <- cleaned_virulence %>% 
  count(vfid, sort = TRUE)

# Preview the top 10 most dominant pathogenic features
head(virulence_counts, 10)

# Plot the distribution of percent identity values 
ggplot(cleaned_virulence, aes(x = pident)) +
  geom_density(fill = "steelblue", alpha = 0.7, color = "darkblue", size = 1) +
  theme_minimal() +
  labs(
    title = "Distribution of Percent Identity for Virulence Factors",
    subtitle = "E. coli Strain N31 (Filtered: Identity ≥ 80%, E-value ≤ 1e-5)",
    x = "Percent Identity (%)",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 11)
  )


# PART 4: Retrieve Gene Nomenclature
# Map R data with the VFDB reference table
# Read the VFDB protein FASTA file
vfdb_fasta <- readLines("data/vfdb/VFDB_setA_pro.fas")

# Filter header lines using useBytes = TRUE for safety
headers <- vfdb_fasta[grep("^>", vfdb_fasta, useBytes = TRUE)]

# Construct a clean gene reference database table
vfdb_ref <- tibble(header = headers) %>%
  mutate(
    vfid = str_extract(header, "VFG[0-9]+"),
    gene_name = str_extract(header, "(?<=\\)).*") %>% str_trim()
  ) %>%
  filter(!is.na(vfid)) %>%
  select(vfid, gene_name)

# Preview the resulting reference mapping database
head(vfdb_ref, 5)


# PART 5: Merge Pathogenic Targets with Clean Gene Nomenclature
# Main inputs:
# - cleaned_virulence: Contains filtered clean pathogenic features from Strain N31
# - vfdb_ref: Dictionary containing translated official gene descriptions from VFDB

# Merge the BLAST results with official VFDB gene names
annotated_virulence <- cleaned_virulence %>%
  inner_join(vfdb_ref, by = "vfid")

# Preview the annotated virulence profile of Strain N31
annotated_virulence %>%
  select(qseqid, vfid, gene_name, pident, evalue) %>%
  head(10)


# PART 6: Prepare Genomic Coordinates (Prokka GFF Processing)
# Define the path to the structural annotation GFF file
gff_path <- "results/annotation/strain_N31.gff"

# Read the GFF file 
# Exclude GFF3 metadata directives (lines prefixed with '##')
gff_data <- read_tsv(
  file = gff_path,
  comment = "##",
  col_names = c("seqid", "source", "type", "start", "end", 
                "score", "strand", "phase", "attributes")
)

# Retain only CDS features and extract the locus identifier from the attributes field
gff_cleaned <- gff_data %>%
  filter(type == "CDS") %>%
  mutate(qseqid = str_extract(attributes, "(?<=ID=)[^;]+")) %>%
  select(qseqid, start, end, strand)

# Preview parsed physical coordinate data
head(gff_cleaned, 5)


# PART 7: Integrate Virulence Identity with Genomic Coordinates
# Combine Gene Name + BLAST Validity Scores + Precise Physical Locus Coordinates

# Merge the annotated virulence data with physical chromosomal positions
final_virulence_map <- annotated_virulence %>%
  inner_join(gff_cleaned, by = "qseqid")

# Preview the comprehensive integrated mapping dataframe
final_virulence_map %>%
  select(qseqid, gene_name, start, end, pident) %>%
  head(5)


# PART 8: Initialize Circular Genomic Mapping (circlize)
# Install and load the circlize engine package
install.packages("circlize")
library(circlize)

# Determine the maximum genomic coordinate to define the chromosome boundary
genome_length <- max(gff_cleaned$end)
genome_length

# Sort features by ascending genomic start position
final_virulence_map_sorted <- final_virulence_map %>%
  arrange(start)

# Dynamic Column Detection
# Dynamically locate the target column containing valid gene nomenclature profiles
gene_column_name <- NA
column_target <- c("g_name", "gene", "target", "subject", "Query_Def", "Subject_Def")

for (column in column_target) {
  if (column %in% colnames(final_virulence_map_sorted)) {
    gene_column_name <- column
    break
  }
}

# Reset circular plot parameters and remove inter-sector gaps
circos.clear()
circos.par(
  "start.degree" = 90,        
  "gap.degree" = 0
)

# Initialise a single circular sector spanning the full genome length
circos.initialize(sectors = "Strain_N31", xlim = c(0, genome_length))

# --- TRACK 1: COORDINATE SCALE (OUTERMOST RING TRACK) ---
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.axis(
    h = "top", 
    labels.cex = 0.65, 
    labels = seq(0, genome_length, by = 100000) / 1000, 
    major.at = seq(0, genome_length, by = 100000),
    labels.facing = "outside"
  )
}, bg.border = NA)

#--- TRACK 2: VIRULENCE FACTORS BLOCK (BORDERED CORE BLOCK TRACK) ---
circos.track(ylim = c(0, 1), bg.col = "#f0f2f5", bg.border = "black")

for(i in 1:nrow(final_virulence_map_sorted)) {
  g_start <- final_virulence_map_sorted$start[i]
  g_end <- final_virulence_map_sorted$end[i]
  g_pident <- final_virulence_map_sorted$pident[i]
  
  block_color <- if(g_pident >= 95) "#d95f02" else "#1f78b4"
  
  circos.rect(
    xleft = g_start, ybottom = 0, xright = g_end, ytop = 1, 
    col = block_color, border = NA, sector.index = "Strain_N31"
  )
}

# --- TRACK INTERMEDIATE: VIRULENCE FACTOR TEXT LABELS (BOUNDED TRACK) ---
# Runs dynamically only if a valid gene nomenclature column is detected
if (!is.na(gene_column_name)) {
  for(i in 1:nrow(final_virulence_map_sorted)) {
    g_start <- final_virulence_map_sorted$start[i]
    g_end <- final_virulence_map_sorted$end[i]
    g_pident <- final_virulence_map_sorted$pident[i]
    
    # Retrieve nomenclature string dynamically using verified column indices
    raw_name <- as.character(final_virulence_map_sorted[i, gene_column_name])
    
    # Text clean-up parsing: extract string within parentheses; default to raw if empty
    short_name <- str_extract(raw_name, "(?<=\\()[^\\)]+")
    if(is.na(short_name) || length(short_name) == 0) short_name <- raw_name
    
    # Render gene labels only for hits with percent identity ≥ 99% to minimise label overlap
    if(!is.na(short_name) & length(short_name) > 0 & g_pident >= 99) {
      circos.text(
        x = (g_start + g_end) / 2, 
        y = 1.7, 
        labels = short_name, 
        sector.index = "Strain_N31",
        facing = "clockwise", 
        niceFacing = TRUE, 
        cex = 0.45, 
        font = 2, 
        col = "darkred"
      )
    }
  }
}

# --- TRACK 3: SKEW/BITSCORE WAVE (CONTINUOUS DENSE WAVE TRACK) ---
min_bit <- min(final_virulence_map_sorted$bitscore, na.rm = TRUE)
max_bit <- max(final_virulence_map_sorted$bitscore, na.rm = TRUE)
mean_bit <- mean(final_virulence_map_sorted$bitscore, na.rm = TRUE)

circos.track(ylim = c(min_bit - 20, max_bit + 20), track.height = 0.18, bg.col = "#f8f9fa", bg.border = "black")

# Construct coordinate vectors for a closed circular area plot
x_wave <- c(0, final_virulence_map_sorted$start, genome_length)
y_wave <- c(mean_bit, final_virulence_map_sorted$bitscore, mean_bit)

circos.lines(
  sector.index = "Strain_N31",
  x = x_wave,
  y = y_wave,
  col = "#2ca02c", 
  type = "l", 
  area = TRUE, 
  baseline = mean_bit
)
  
# Append the map legend
text(0, 0, "Refined Genomic Mapping\nE. coli Strain N31\n(Complete Circular View)", cex = 0.85, font = 2)
legend(
  "topright", 
  legend = c("Identity ≥ 95%", "Identity < 95%", "Bitscore Waves"), 
  fill = c("#d95f02", "#1f78b4", "#2ca02c"),
  border = NA,
  bty = "n",         
  cex = 0.8          
)


# PART 9: Profiling Top 20 Virulence Factors
library(ggplot2)
library(dplyr)
library(stringr)

# Extract gene nomenclature from the annotated virulence map
top_data_prep <- final_virulence_map %>%
  mutate(raw_gene = as.character(gene_name)) %>%
  
# Clean up gene text identifiers inside parentheses (e.g., parsing 'stx1' out of complex rows)
mutate(clean_gene_name = str_extract(raw_gene, "(?<=\\()[^\\)]+")) %>%
mutate(clean_gene_name = ifelse(is.na(clean_gene_name) | clean_gene_name == "", raw_gene, clean_gene_name))

# Filter the top 20 virulence determinants
top_virulence <- top_data_prep %>%
  arrange(desc(pident), desc(bitscore)) %>%  
  head(20)

# Visualise the top 20 virulence factors using ggplot2
ggplot(top_virulence, aes(x = reorder(clean_gene_name, pident), y = pident, fill = pident)) +
  geom_col(width = 0.65, color = "#2b2b2b", linewidth = 0.2) + 
  coord_flip(ylim = c(min(top_virulence$pident) - 0.5, 100)) + 
  
scale_fill_gradientn(
  colors = c("#3182bd", "#e6550d", "#b30000"),
  labels = function(x) paste0(x, "%")
) +
  
theme_minimal(base_family = "sans") +
  labs(
  title = "Top 20 Highly Identical Virulence Factors",
  subtitle = "E. coli Strain N31 (PacBio Long-Reads Analysis)",
  x = "Virulence Gene Factor",
  y = "Percent Identity (%)",
  fill = "Identity"
) +
  
theme(
  plot.title = element_text(face = "bold", size = 14, color = "#1a1a1a", margin = margin(b = 5)),
  plot.subtitle = element_text(size = 10, face = "italic", color = "#555555", margin = margin(b = 15)),
  axis.text.y = element_text(size = 9.5, face = "bold.italic", color = "#2c3e50"), 
  axis.text.x = element_text(size = 9, color = "#2c3e50"),
  axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 10)),
  axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 10)),
  panel.grid.major.y = element_blank(), 
  panel.grid.minor = element_blank(),
  legend.position = "right",
  legend.title = element_text(size = 9, face = "bold")
)