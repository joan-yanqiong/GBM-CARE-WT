#################################
## Title: Extended Data Figure 7d in Nomura et al
## Date: 2025.02.13
## Author: Kevin Johnson
## Description: Produce heatmap for tumor samples that were profiled by both bulk RNA sequencing and single nucleus RNA seq. Note samples were taken from different sections of the same tumor.
#################################

# Necessary packages
library(tidyverse) # v1.3.1
# Custom minimalist plotting theme
source("R/plot_theme.R")

# Data including bulk RNA TCGA subtype assignment (based on highest enrichment and p-value) and snRNA baseline profile assignment.
care_scp_rna <- readRDS("data/care_snrna_bp_bulk_rna_subtype.RDS")

## Test whether state controlled single nucleus Baseline Profile is associated with bulk RNA TCGA subtype.
# Perform analysis both with and without "Mixed" class to see whether it has an impact on results.
care_scp_rna_no_mixed <- care_scp_rna %>% 
  filter(snrna_bp!="Mixed")

# Both Fisher exact p-values are  p = 0.02. Excluding Mixed does not impact the resulting association.
fisher.test(table(care_scp_rna$snrna_bp, care_scp_rna$bulk_tcga_subtype))
fisher.test(table(care_scp_rna_no_mixed$snrna_bp, care_scp_rna_no_mixed$bulk_tcga_subtype))

# Format data for plotting.
bp_bulk_plot <- care_scp_rna %>% 
  group_by(snrna_bp, bulk_tcga_subtype) %>% 
  summarise(counts = n()) %>% 
  dplyr::select(snrna_bp, bulk_tcga_subtype, counts)

# Set order for plotting.
bp_bulk_plot$snrna_bp <- factor(bp_bulk_plot$snrna_bp, levels = c("BP-ECM", "BP-Neuronal", "BP-Glial", "Mixed"))

pdf("figures/nomura_edf7d.pdf", width = 6, height = 4, useDingbats = FALSE)
ggplot(bp_bulk_plot, aes(snrna_bp, bulk_tcga_subtype)) +
  geom_tile(aes(fill = counts)) +  
  geom_text(aes(label = counts)) +
  scale_fill_gradient(low = "white", high = "red") +
  plot_theme +
  theme(legend.position = "top", 
        axis.text.x = element_text(angle=90, hjust=1)) +
  labs(y = "Bulk RNA TCGA subtype", x = "snRNA state controlled baseline profile")
dev.off()

### END ###