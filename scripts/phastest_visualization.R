# ==============================================================================
# PROPHAGE REGION VISUALIZATION PIPELINE (Escherichia coli Strain N31)
# Description: This script generates both linear (ggplot2) and circular (circlize) 
#              maps representing the distribution, completeness, and functional 
#              annotation of prophage regions predicted by PHASTEST.
# ==============================================================================

# 1. ENVIRONMENT SETUP & DEPENDENCIES
# Uncomment the line below if you need to install the ggrepel package first:
# install.packages("ggrepel")

library(ggplot2)
library(dplyr)
library(ggrepel)
library(circlize)

# 2. ESSENTIAL DATASET (PHASTEST Predictions)
phastest_data <- data.frame(
  Region = 1:12,
  Length_Kb = c(31.9, 45.2, 35.8, 32.6, 42.5, 42.6, 28.9, 9.8, 9.8, 7.8, 7.0, 7.6),
  Completeness = c("intact", "intact", "questionable", "questionable", "intact", "intact", 
                   "intact", "intact", "intact", "questionable", "incomplete", "incomplete"),
  Score = c(150, 150, 80, 70, 130, 150, 150, 120, 120, 80, 40, 50),
  Start = c(443643, 700630, 1018182, 1453961, 2447567, 2535491, 3515673, 3599613, 3609461, 4564053, 5001927, 5171645),
  End = c(475583, 745884, 1054032, 1486654, 2490133, 2578119, 3544667, 3609458, 3619295, 4571914, 5008981, 5179244),
  Most_Common_Phage = c("Klebsi_4LV2017", "Burkho_BcepMu", "Entero_SfI", "Entero_mEp460", 
                        "Entero_Sf101", "Entero_mEp460", "Shigel_SfII", "Entero_lambda", 
                        "Entero_HK630", "Escher_500465_1", "Entero_YYZ_2008", "Entero_lambda")
)

genome_length <- 5179245 

# Factorize completeness for consistent color mapping scale
phastest_data$Completeness <- factor(phastest_data$Completeness, 
                                     levels = c("intact", "questionable", "incomplete"))

# Calculate midpoint positions for label placement
phastest_data$Mid <- (phastest_data$Start + phastest_data$End) / 2
phastest_data$Label_Full <- paste0("R", phastest_data$Region, ": ", phastest_data$Most_Common_Phage)


# 3. LINEAR PROPHAGE MAP (ggplot2)
linear_plot <- ggplot(phastest_data) +
  # Genome backbone
  geom_hline(yintercept = 0, color = "#7f8c8d", linewidth = 2) +
  
  # Prophage region bars
  geom_rect(aes(xmin = Start, xmax = End, ymin = -0.15, ymax = 0.15, fill = Completeness), 
            color = "#2c3e50", linewidth = 0.3, alpha = 0.9) +
  
  # Top labels: Region ID & Most Common Phage Hit
  geom_text_repel(
    aes(x = Mid, y = 0.2, label = paste0("R", Region, ": ", Most_Common_Phage)),
    direction = "y",           
    nudge_y = 0.4,             
    force = 2,                 
    box.padding = 0.3, 
    segment.color = "grey50",  
    segment.size = 0.4,        
    size = 3.2, 
    fontface = "bold"
  ) +
  
  # Bottom labels: Region Size & PHASTEST Score
  geom_text_repel(
    aes(x = Mid, y = -0.2, label = paste0(Length_Kb, " Kb\n(Score: ", Score, ")")),
    direction = "y",
    nudge_y = -0.4,            
    force = 1.5,
    box.padding = 0.2,
    segment.color = "grey50",
    segment.size = 0.4,
    size = 2.6, 
    lineheight = 0.85          
  ) +
  
  # Color customization
  scale_fill_manual(
    values = c("intact" = "#2ca02c", "questionable" = "#ff7f0e", "incomplete" = "#d62728"),
    labels = c("Intact (Score > 90)", "Questionable (Score 70-90)", "Incomplete (Score < 70)")
  ) +
  
  # X-axis configuration (Converting to Mb)
  scale_x_continuous(
    limits = c(0, genome_length),
    labels = function(x) paste0(x / 1e6, " Mb"),
    breaks = seq(0, genome_length, by = 1e6)
  ) +
  
  ylim(-1.5, 1.5) +
  
  labs(
    title = "Distribution of Identified Prophage Regions across the Genome",
    subtitle = "Escherichia coli Strain N31 (PHASTEST Functional Annotation)",
    x = "Genome Position (Mb)",
    y = "",
    fill = "Completeness Status"
  ) +
  
  theme_minimal(base_family = "sans") + 
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(color = "#eef0f2", linewidth = 0.5), 
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14, color = "#2c3e50"),
    plot.subtitle = element_text(face = "italic", size = 10, color = "#7f8c8d", margin = margin(b = 10))
  )

# Display the linear plot
print(linear_plot)


# 4. CIRCULAR GENOME MAP (circlize)
# Reset layout parameters
circos.clear()
circos.par(
  start.degree = 90, 
  gap.degree = 0, 
  track.margin = c(0.01, 0.01),
  points.overflow.warning = FALSE
)

# Initialize circular layout with a single chromosome track
circos.initialize(factors = "chr1", xlim = c(0, genome_length))

# --- TRACK 1: STAGGERED LABELS (Prevents text overlap in dense clusters) ---
circos.track(factors = "chr1", ylim = c(0, 1), track.height = 0.35, bg.border = NA,
             panel.fun = function(x, y) {
               for(i in 1:nrow(phastest_data)) {
                 
                 # Interlocking staggered formula for Alternating Heights
                 if (i %% 2 == 0) {
                   y_text <- 0.05       # Text placement for lower row (closer to ring)
                   y_line_end <- 0.03   # Short leader line
                 } else {
                   y_text <- 0.48       # Text placement for upper row (pushed outwards)
                   y_line_end <- 0.46   # Long leader line
                 }
                 
                 # Draw labels
                 circos.text(x = phastest_data$Mid[i], 
                             y = y_text, 
                             labels = phastest_data$Label_Full[i], 
                             facing = "clockwise", 
                             niceFacing = TRUE, 
                             cex = 0.65, 
                             font = 2, 
                             adj = c(0, 0.5), 
                             col = "#2c3e50")
                 
                 # Draw leader lines
                 circos.lines(x = c(phastest_data$Mid[i], phastest_data$Mid[i]), 
                              y = c(0, y_line_end), 
                              col = "#7f8c8d", lwd = 1)
               }
             })

# --- TRACK 2: CHROMOSOME SCALE AXIS (Position markers in Mb) ---
circos.track(factors = "chr1", ylim = c(0, 1), track.height = 0.06, bg.border = NA,
             panel.fun = function(x, y) {
               circos.axis(h = "bottom", labels.cex = 0.65, labels.font = 2,
                           major.at = seq(0, genome_length, by = 1e6),
                           labels = paste0(seq(0, genome_length, by = 1e6)/1e6, " Mb"),
                           col = "#34495e", labels.col = "#2c3e50")
             })

# --- TRACK 3: MAIN GENOME RING & PROPHAGE BARS ---
circos.track(factors = "chr1", ylim = c(0, 1), track.height = 0.15, 
             bg.col = "#f8f9fa", bg.border = "#2c3e50")

# Map colored rectangles inside the main ring base on completeness status
for(i in 1:nrow(phastest_data)) {
  color_choice <- switch(as.character(phastest_data$Completeness[i]),
                         "intact" = "#2ca02c",
                         "questionable" = "#ff7f0e",
                         "incomplete" = "#d62728")
  circos.rect(phastest_data$Start[i], 0, phastest_data$End[i], 1, 
              sector.index = "chr1", col = color_choice, border = "#1e272e")
}

# Add Main Diagram Title
title(main = "Prophage Distribution Map (E. coli Strain N31)", cex.main = 1.2, font.main = 2)

# --- MAP LEGEND ---
legend(x = "bottomright", 
       legend = c("Intact (Score > 90)", "Questionable (Score 70-90)", "Incomplete (Score < 70)"), 
       fill = c("#2ca02c", "#ff7f0e", "#d62728"), 
       bty = "n", 
       cex = 0.75,           
       pt.cex = 1.2,        
       title = "Completeness Status", 
       title.font = 2)
