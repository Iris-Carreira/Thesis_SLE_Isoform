library(IsoformSwitchAnalyzeR)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Set variables for the analysis
working_dir <- "C:/Users/paulo/OneDrive - FCT NOVA/Desktop/postdoc/project_LUPUS/public_datasets/IsoformSwitchAnalyzeR/GSE118254/"

# Data Import
setwd(working_dir)

# Ensure the object name matches throughout the script
switch_list_stats <- readRDS("IsoformSwitchResults_SM.rds")

# Identify top hits for labeling based on user-defined thresholds
top_hits <- SwitchAnalysisResults_PD1POS_Real$isoformFeatures %>%
  filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1) %>%
  group_by(condition_2) %>% # Grouping ensures labels appear for each facet
  #slice_min(order_by = isoform_switch_q_value, n = 20, with_ties = FALSE)
  slice_max(order_by = abs(dIF), n = 20, with_ties = FALSE)

# --- 1. Visualize Volcano Plot: Isoform Switch Significance ---

ggplot(data = SwitchAnalysisResults_PD1POS_Real$isoformFeatures, aes(x = dIF, y = -log10(isoform_switch_q_value))) +
  geom_point(
    aes(fill = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05), # Use 'fill' for the center
    shape = 21,                # Circle with a border
    color = "white",           # The border color
    size = 3,                # Slightly larger to show off the stroke
    stroke = 0.4,              # The thickness of the white line
    alpha = 0.5                # Transparency for the whole marker
  ) +
  geom_text_repel(
    data = top_hits, 
    aes(label = gene_name), 
    size = 3, # Slightly smaller text usually helps with clutter
    box.padding = 0.5, 
    max.overlaps = Inf, # Allows more labels to show if they are tight
    segment.color = 'grey'
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 'dashed') +
  facet_wrap(~ condition_2) +
  scale_fill_manual('Significant Switch', values = c('darkgrey', 'salmon')) + 
  labs(title = "Isoform Switch Volcano Plot", x = 'dIF', y = '-Log10(Q-Value)') +
  theme_bw()

ggsave("Volcano.Isoforms.PD1POS.png", width = 10, height = 8, dpi = 300)

# --- 2. Visualize Volcano Plot: Isoform Switch vs. Gene Change ---

# plotting
ggplot(data = SwitchAnalysisResults_PD1POS_Real$isoformFeatures, aes(x = gene_log2_fold_change, y = dIF)) +
  geom_point(
    aes(fill = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05), # Use 'fill' for the center
    shape = 21,                # Circle with a border
    color = "white",           # The border color
    size = 3,                  # Slightly larger to show off the stroke
    stroke = 0.4,              # The thickness of the white line
    alpha = 0.5                # Transparency for the whole marker
  ) +
  geom_text_repel(
    data = top_hits, 
    aes(label = gene_name), 
    size = 3, # Slightly smaller text usually helps with clutter
    box.padding = 0.5, 
    max.overlaps = Inf, # Allows more labels to show if they are tight
    segment.color = 'grey'
  ) +
  facet_wrap(~ condition_2) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_fill_manual('Significant Switch', values = c('darkgrey', 'salmon')) + 
  labs(
    title = "Isoform Switch vs Gene Expression Change",
    x = 'Gene Log2 Fold Change', 
    y = 'dIF (Delta Isoform Fraction)'
  ) +
  theme_bw()

ggsave("Volcano.Isoforms.vs.GeneExpression.PD1POS.png", width = 10, height = 8, dpi = 300)

# Overview of a given gene
goi <- 'XIST'

# Use paste0 to insert the variable into the path
png(file = paste0('switchPlot_', goi, '_Tfh.png'), height = 500, width = 800, res = 100)
switchPlot(switch_list_stats, gene = goi, plotTopology = FALSE)
dev.off()


switch_list_stats$switchConsequence


