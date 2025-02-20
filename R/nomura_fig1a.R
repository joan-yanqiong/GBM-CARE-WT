#################################
## Title: Figure 1 panel a in Nomura et al - genetic alteration frequency
## Date: 2025.02.15
## Author: Kevin Johnson
## Description: Produce a cohort-level breakdown for the relative frequency of different genetic alteration from accompanying bulk WES/WGS data
#################################

# Necessary packages
library(tidyverse) # v1.3.1
source("R/plot_theme.R")

# Patient-level longitudinal information. tumor_a and tumor_b represent the first two time points for a patient/case.
# To determine arm gain/loss requires that a CNA segment be consistently gained or lost. Samples with disrupted segments where some proportions are gain while others are lost are assigned an NA.
arm_heatmap <- readRDS("data/genetic_patient_chr710_arm_status.RDS")
# Selected gene-level copy number alterations
cnv_heatmap <- readRDS("data/genetic_patient_gene_cna.RDS")
# Selected gene-level driver gene mutations
mut_heatmap <- readRDS("data/genetic_patient_gene_mutation.RDS")
# Cohort-level summary for both mutation and copy number alterations by time point.
snvgdata <- readRDS("data/genetic_cohort_gene_mutation.RDS")
cnvgdata <- readRDS("data/genetic_cohort_gene_cna.RDS")

# Select the features and order for the alterations of interest. Selected features based on PMID: 31748746
snv_gene_order <- c("TERT", "PTEN", "EGFR", "TP53", "NF1", "RB1", "PIK3R1")
cnv_gene_order <- c("CDKN2A", "EGFR", "CDK4", "PDGFRA", "MDM2")
arm_order <- rev(c("7q","10q"))

# Organize the data as a summary breakdown across time points.
mut_values <- snvgdata %>%
  gather(key = "type", value = "value", shared, private_a, private_b) %>%
  mutate(type = factor(type,
                       levels = c("private_b", "private_a", "shared"),
                       labels = c("T2", "T1", "T1&T2")),
         # Total number of patients
         pct_value = value/n_distinct(mut_heatmap$patient),
         gene_symbol = factor(gene_symbol,
                              levels = c("TERT", "PTEN", "EGFR", "TP53", "NF1", "RB1", "PIK3R1"),
                              labels = c("TERTp", "PTEN", "EGFR","TP53", "NF1", "RB1", "PIK3R1"))) 

cna_values <- cnvgdata %>%
  gather(key = "type", value = "value", shared, private_a, private_b) %>%
  mutate(type = factor(type,
                       levels = c("private_b", "private_a", "shared"),
                       labels = c("T2", "T1", "T1&T2"))) %>% 
  mutate(gene_symbol = factor(gene_symbol,
                              levels = c("CDKN2A", "EGFR", "MDM2","CDK4", "PDGFRA"),
                              labels = c("CDKN2A del", "EGFR amp",  "MDM2 amp", "CDK4 amp", "PDGFRA amp")),
         pct_value = value/n_distinct(cnv_heatmap$patient)) 

# For chromosome arms, I am using 7q as proxy for Chr7 and 10q as proxy for Chr10 because there is too much information loss with the chr7 and chr10 query. 
# That is, there is a missingness of 10 cases due to disrupted segments. In that analysis, 27 out of 36 (75%) had Chr7+/10- co-alterations (10 NA values).
arm_heatmap_extnd <- arm_heatmap %>% 
  filter(arm%in%c("7q","10q")) %>% 
  mutate(segment_qual = ifelse(!is.na(cnv_state), "Contiguous", "Disrupted"),
         cnv_state = ifelse(cnv_state == "neut", NA, cnv_state),
         arm = factor(arm, levels = arm_order))

# Summarizing the contiguous calls. That is, those samples where a more definitive copy number call could be made.
# Note: Our approach doesn't directly take into account tumor purity, which may reduce the number of called chr7/chr10 alterations.
arm_continuous_measurements <- arm_heatmap_extnd %>% 
  group_by(arm, segment_qual) %>% 
  summarise(countinuous_counts = n()) %>% 
  filter(segment_qual=="Contiguous")

arm_processed_values <- arm_heatmap_extnd %>%
  # Restricting to those events that are chr7q amplifications and chr10q deletions.
  mutate(gene_symbol = paste(arm, cnv_state, sep="_")) %>% 
  filter(gene_symbol%in%c("7q_amp", "10q_del")) %>% 
  mutate(type = factor(cnv_change,
                       levels = c("R", "P", "S"),
                       labels = c("T2", "T1", "T1&T2"))) %>%
  mutate(gene_symbol = factor(gene_symbol,
                              levels = c("10q_del", "7q_amp"),
                              labels = c("Chr10 loss", "Chr7 gain"))) %>% 
  dplyr::select(tumor_pair_barcode:tumor_barcode_b, gene_symbol, type) %>% 
  group_by(gene_symbol, type) %>% 
  summarise(value = n()) 

# Since I am tabulating by the long arm of each chromsome and there are differences in the number of disrupted segments. 
# An arm can't have a gain/loss is the segment is disrupted. Note there are more T1-only values, likely due to purity decreased (i.e., decreased malignant cells, increased oligodendrocytes at recurrence).
arm_values <- arm_processed_values %>% 
  mutate(pct_value = case_when(
    gene_symbol == "Chr10 loss" ~ value / arm_continuous_measurements$countinuous_counts[arm_continuous_measurements$arm=="10q"],
    gene_symbol == "Chr7 gain" ~ value / arm_continuous_measurements$countinuous_counts[arm_continuous_measurements$arm=="7q"],
    TRUE ~ NA # Keep the original value if gene_symbol doesn't match either "Chr10" or "Chr7"
  ))

# Format each variant type to be merged altogether.
mut_merge <- mut_values %>% 
  dplyr::select(gene_symbol, type, pct_value) %>% 
  mutate(alt = "SNV")
cna_merge <- cna_values %>% 
  dplyr::select(gene_symbol, type, pct_value) %>% 
  mutate(alt = "CNA") 
arm_merge <- arm_values %>% 
  dplyr::select(gene_symbol, type, pct_value) %>% 
  mutate(alt = "Arm-level")

all_summary <- bind_rows(mut_merge, cna_merge, arm_merge)
all_summary$gene_symbol <- factor(all_summary$gene_symbol, levels=c("Chr10 loss", "Chr7 gain", "CDKN2A del", "EGFR amp", "MDM2 amp","CDK4 amp",  "PDGFRA amp", "TERTp", "PTEN", "EGFR", "TP53", "NF1", "RB1", "PIK3R1"))

all_summary_df <- all_summary %>%
  group_by(gene_symbol, alt) %>%
  summarise(gene_alt_pct_value = sum(pct_value)) %>% 
  ungroup()

pdf("figures/nomura_fig1b_genetic_alteration_overview.pdf", width = 6, height = 4, useDingbats = FALSE)
all_summary_df %>% ggplot(aes(x = gene_symbol, y = gene_alt_pct_value*100)) +
  geom_bar(position="stack", stat="identity", fill = "#377EB8") +
  labs(y = "Percent of patients\nwith alteration (%)", x = "Chromosome/gene feature", fill="Variant Evolution") +
  guides(fill="none") +
  ylim(0, 100) + 
  facet_grid(.~alt, scales = "free", space = "free") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

### END ###