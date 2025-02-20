#################################
## Title: Spitzer Figure 1d - genetic landscape across time points
## Date: 2025.02.18
## Author: Kevin Johnson
## Description: Generate donut chart for average mutation fraction and driver event landscape panel across patients
#################################

# Necessary packages
library(tidyverse) # v1.3.1
source("R/plot_theme.R")

# Inspecting mutations where there is coverage of at least 14x to determine relative longitudinal mutation fraction.
mut_freq_prop_case <- readRDS("data/genetic_patient_mutation_fraction.RDS")

# To determine arm gain/loss requires that a CNA segment be consistently gained or lost. Samples with disrupted segments where some proportions are gain while others are lost are assigned an NA.
arm_heatmap <- readRDS("data/genetic_patient_chr710_arm_status.RDS")
# Selected gene-level copy number alterations
cnv_heatmap <- readRDS("data/genetic_patient_gene_cna.RDS")
# Selected gene-level driver gene mutations
mut_heatmap <- readRDS("data/genetic_patient_gene_mutation.RDS")

# Cohort-level summary for both mutation and copy number alterations by time point.
snvgdata <- readRDS("data/genetic_cohort_gene_mutation.RDS")
cnvgdata <- readRDS("data/genetic_cohort_gene_cna.RDS")

### ### ### ### ### ###
### Top - Generate mutation proportion donut chart
### ### ### ### ### ###
res_mut_prop <- mut_freq_prop_case %>% 
  mutate(mutation_type = recode(mutation_type, `P` = "T1",
                                `R` = "T2",
                                `S` = "T1&T2")) %>% 
  group_by(mutation_type) %>% 
  summarise(avg_prop = mean(mutation_percent),
            med_prop = median(mutation_percent))

data <- data.frame(
  category=c("T2-only", "T1+T2", "T1-only"),
  avg_mean=c(res_mut_prop$avg_prop[res_mut_prop$mutation_type=="T2"],  res_mut_prop$avg_prop[res_mut_prop$mutation_type=="T1&T2"], res_mut_prop$avg_prop[res_mut_prop$mutation_type=="T1"])
)

data$avg_mean <- round(data$avg_mean, digits = 3)
data$category <- factor(data$category, levels=c("T2-only", "T1+T2", "T1-only"))

# Compute fractions
data$fraction <- data$avg_mean / sum(data$avg_mean)

data$ymax <- cumsum(data$fraction)

data$ymin <- c(0, head(data$ymax, n=-1))

data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a descriptive label
data$label <- paste0(data$category, "\n", data$avg_mean*100, " %")

pdf("figures/spitzer_fig1d_mutation_burden_donut.pdf", width = 1.5, height = 1.5)
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_manual(values=c("T1+T2" = "#CA932F",
                             "T2-only" = "#2FB3CA",
                             "T1-only" = "#CA2F66")) +
  coord_polar(theta = "y", start = 0)+
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")
dev.off()


### ### ### ### ### ### 
### Fig. 1d - bottom panel
### ### ### ### ### ###
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

# Sample fraction with de novo HM (2/(54 NB+tonly)), treatment assoc. HM (4/(52 NB+tonly - 2 de novo)), and acquired small deletion ((10/45 NB - 1 de novo))
phenotype_merge <- data.frame(
  gene_symbol = c(rep("De novo HM", 3), rep("Treatment assoc. HM", 3), rep("Small deletion", 3)),
  type = c("T2", "T1", "T1&T2","T2", "T1", "T1&T2", "T2", "T1", "T1&T2"),
  pct_value = c(0, 0, .04, .08, 0, 0, .22, 0, 0)
)
phenotype_merge <- phenotype_merge %>% mutate(alt = "Phenotype")
phenotype_merge$type <- as.factor(phenotype_merge$type)
all_summary <- bind_rows(mut_merge, cna_merge, arm_merge, phenotype_merge)

all_summary$gene_symbol <- factor(all_summary$gene_symbol, levels=c("Chr10 loss", "Chr7 gain", "CDKN2A del", "EGFR amp", "MDM2 amp","CDK4 amp",  "PDGFRA amp", "TERTp", "PTEN", "EGFR", "TP53", "NF1", "RB1", "PIK3R1", "De novo HM", "Treatment assoc. HM", "Small deletion"))
all_summary$alt <- factor(all_summary$alt, levels=c("Arm-level", "CNA", "SNV", "Phenotype"))

pdf("figures/spitzer_fig1d_genetics_overview.pdf", width = 7, height = 5, useDingbats = FALSE, bg="transparent")
all_summary %>% ggplot(aes(x = gene_symbol, y = pct_value*100, fill = type)) +
  geom_bar(position="stack", stat="identity") +
  labs(y = "Percent of cases with alteration (%)", x = "", fill="Variant Evolution") +
  guides(fill="none") +
  scale_fill_manual(values=c("T2" = "#2FB3CA", "T1" ="#CA2F66", "T1&T2"="#CA932F")) +
  ylim(0, 100) + 
  facet_grid(.~alt, scales = "free", space = "free") +
  plot_theme +
  theme(strip.background = element_rect(fill = "white", color = "white")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Note: combine top and bottom panel in illustrator

### END ###