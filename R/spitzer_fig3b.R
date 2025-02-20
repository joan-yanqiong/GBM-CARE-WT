#################################
## Title: Figure 3b in Spitzer et al - mutation fraction panel for patients with 3 time points
## Date: 2025.02.15
## Author: Kevin Johnson
## Description: Produce a stacked mutation fraction bar plot for 2 patients with 3 time points
#################################

# Necessary packages
library(tidyverse) # v1.3.1
source("R/plot_theme.R")


# File containing information to link snrna and bulk dna identifiers.
linker_file_samples <- read.table("data/care_wt_synapse_dna_snrna_id_linker.txt", sep="\t", header = TRUE)
linker_file_cases <- linker_file_samples %>% 
  mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  dplyr::select(Patient, case_barcode) %>% 
  distinct()
  
# Breaks down the mutation count for longitudinal comparisons. That is, which mutations are detected at one time point and which ones are shared across multiple time points.
tumor_mut_comparison <- read.delim("data/tumor_mut_comparison.txt", header = TRUE, sep = "\t")
tumor_mut_comparison <- tumor_mut_comparison %>% 
  left_join(linker_file_cases, by="case_barcode")


# Create the data.frame for the CARE patients P20 and P32
mut_freq_prop_case = tumor_mut_comparison %>% 
  filter(Patient%in%c("P20", "P32")) %>% 
  mutate(P = (count_a-intersection_ab)/union_ab,
         R =  (count_b-intersection_ab)/union_ab,
         S = intersection_ab/union_ab) %>% 
  select(P, R, S, tumor_pair_barcode, union_ab) %>% 
  gather(mutation_type, mutation_percent, c(P, R, S), -tumor_pair_barcode, -union_ab) %>%
  mutate(mutation_type = factor(mutation_type, levels = c("R", "P", "S"))) %>% 
  mutate(tumor_comparison_a = substr(tumor_pair_barcode, 14, 15),
         tumor_comparison_b = substr(tumor_pair_barcode, 20, 21),
         tumor_comparison = paste0(tumor_comparison_a, "-", tumor_comparison_b),
         case_barcode = substr(tumor_pair_barcode, 9, 12))
mut_freq_prop_case$tumor_comparison <- factor(mut_freq_prop_case$tumor_comparison, levels = c("TP-R1", "R1-R2", "TP-R2"))

# Combine the mutation fraction panel with clinical trajectory and relative cell abundance and baseline profile score in Illustrator.
pdf("figures/spitzer_fig3b_mutation_fraction_panel.pdf",width=6, height=3.5, useDingbats = FALSE)
ggplot(mut_freq_prop_case %>%
         mutate(t_comparison = recode(tumor_comparison, `TP-R1` = "T1-T2",
                                      `R1-R2` = "T2-T3",
                                      `"TP-R2` = "T1-T3"),
                patient = recode(case_barcode, `TK01` = "P20",
                                 `TK13` = "P32")) %>% 
         filter(tumor_comparison%in%c("TP-R1", "R1-R2")), aes(x = t_comparison, y = mutation_percent, fill=mutation_type)) +
  geom_bar(stat="identity") +
  labs(y = "% mutation") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F"),
                    labels=c("R" = "Private-Late", "P" ="Private-Early", "S"="Shared")) +
  facet_grid(.~patient, scales="free", space="free") +
  labs(y="Mutation fraction", x="", fill="Mutation Status") +
  plot_theme
dev.off()
