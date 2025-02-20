#################################
## Title: Figure 5 panel a in Nomura et al - genetic alterations and malignant cell state abundance
## Date: 2025.02.20
## Author: Kevin Johnson
## Description: Produce a box plot of selected malignant cell state abundance across groups that either have or do not have a common GBM genetic alteration
#################################

# Necessary packages
library(tidyverse) # v1.3.1
library(ggpubr) # v0.6.0
library(EnvStats) # v2.8
source("R/plot_theme.R")

# Read in processed genetic data and metadata
# File that maps snRNA IDs to aliquot_barcode (identifiers for DNA sequencing).
linker_file_samples <- read.table("data/care_wt_synapse_dna_snrna_id_linker.txt", sep="\t", header = TRUE)

# Malignant cell annotation used throughout the manuscript.
malignant <- readRDS("data/malignant_meta_data_2025_01_08.RDS")

# Gene-level copy number estimation. "hlvl_call" = High-level copy number call. Possible values include: -2 = likely homozygous (high-level) deletion, -1 = likely deletion, 0 = likely copy neutral, 1 = likely amplification, 2 = likely high-level amplification.
gene_copy_number_care <- readRDS("data/genetic_select_gene_level_cna.RDS")

# Chromosome arm call for loss (-1), neutral (0), gain (1) or not able to determine (NA). NA may result from disrupted segments where one part of the arm may be amplified while another deleted.
gatk_cnv_by_arm_care <- readRDS("data/genetic_sample_level_chr_arm_cna.RDS")

# Mutation table to determine whether there was a mutation detected (1) or not (0) in a particular sample. This is restricted to those tumor samples with matching blood.
care_snv_df <- readRDS("data/driver_snv_status_care_20250220.RDS") 

### First, determine the relative abundance of malignant states in each of these tumors.
malignant$Patient <-  sapply(strsplit(malignant$ID, "T"), "[[", 1)
malignant$Timepoint <-  paste0("T", sapply(strsplit(malignant$ID, "T"), "[[", 2))

# Examine malignant frequency. Tabulate cell cycle first.
malignant_freq_cc <- malignant %>% 
  group_by(ID, isCC) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>%
  complete(ID, isCC,
           fill = list(counts = 0, freq = 0)) %>% 
  filter(isCC==TRUE) %>% 
  mutate(State = "Cycling",
         Population = "cell_cycle") %>% 
  dplyr::select(ID, State, counts, freq, Population)

# Combine malignant cell proportions + cycling percentage.
all_cell_freq = malignant %>%
  group_by(ID, State) %>%
  mutate(State = paste0(State, "-like")) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(ID, State,
           fill = list(counts = 0, freq = 0)) %>%
  mutate(Population = "malignant") %>% 
  # Append cycling percentage
  bind_rows(malignant_freq_cc) 

# Remove any sample with less than 50 tumor cells
malignant_cell_number <- malignant %>% 
  group_by(ID) %>% 
  summarise(counts = n()) 

# Set threshold for samples that should be kept in an analysis at 50 malignant cells.
malignant_cell_number_keep <- malignant_cell_number %>% 
  filter(counts > 49)

# Creating a data frame with "high-quality" malignant cells. This is to reduce noise associated with samples that only have a few malignant cells.
all_malignant_freq_hq <- all_cell_freq %>% 
  filter(ID%in%malignant_cell_number_keep$ID) %>% 
  mutate(Patient = sapply(strsplit(ID, "T"), "[[", 1)) %>% 
  mutate(State = recode(State, `Hypoxia-like` = "Hypoxia",
                        `Neuron-like` = "NEU-like",
                        `Stress-like` = "Stress"))


### Reshape some of the genetic tables.
# Restrict to genes with recognized GBM importance and set a cnv_status variable to assign a 0 or 1 depending on whether the sample has the event of interest (e.g., EGFR amplification = 1).
gene_copy_number_care_select <- gene_copy_number_care %>% 
  filter(gene_symbol%in%c("EGFR", "CDK4", "PDGFRA", "MDM2", "CDKN2A")) %>% 
  mutate(cnv_status = ifelse((hlvl_call==2 & gene_symbol%in%c("EGFR", "CDK4", "PDGFRA", "MDM2")), 1, 
                             # Note: Due to reduced purity at recurrence and reduced sensitivity to detect -2 vs -1 in bulk BOTH deletions events (Hom./Het.) are grouped together for this analysis.
                             ifelse((hlvl_call%in%c(-1, -2) & gene_symbol%in%c("CDKN2A")), 1, 0)))

# Setting up an analysis for chromosome 7 gain and 10 loss.
gatk_cnv_by_arm_care_select <- gatk_cnv_by_arm_care %>% 
  filter(arm%in%c("7p", "10p")) %>% 
  # Restricting to a specific arm-level event because there are some non-full chromosome changes.
  mutate(cnv_status = ifelse((arm_call==1 & arm=="7p"), 1, ifelse((arm_call==-1 & arm=="10p"), 1, 0))) %>% 
  # Reformatting to merge with gene-level analyses.
  dplyr::select(ID, Patient, Timepoint, aliquot_barcode, gene_symbol=arm, hlvl_call=arm_call, cnv_status)

# There are 108 samples with both sufficient malignant cells and bulk DNA copy number calls
cnv_state_long_hq <- gene_copy_number_care_select %>% 
  bind_rows(gatk_cnv_by_arm_care_select) %>% 
  inner_join(all_malignant_freq_hq, by=c("ID", "Patient"))

# SNVs 
care_snv_long <- care_snv_df %>% 
  dplyr::select(ID, TERT:PIK3CA) %>% 
  pivot_longer(cols= c(TERT:PIK3CA),
               names_to = "gene_symbol",
               values_to = "snv") 

# There are 87 samples with both SNV calls and sufficient malignant cell number.
snv_state_long_hq <- all_malignant_freq_hq %>% 
  inner_join(care_snv_long, by="ID")


### ### ### ### ### ### ### ###
### Statistics and visualization
### ### ### ### ### ### ### ###
# SNV - filter to mutations significantly associated with cell state abundance. See Neftel et al. and `nomura_edf10a_c.R` for identification of these genes.
care_all_snv_selected <- snv_state_long_hq %>% 
  filter((gene_symbol=="TP53" & State=="Cycling") | 
           (gene_symbol=="TERT" & State=="GPC-like") | 
           (gene_symbol=="TP53" & State=="NEU-like") |
           (gene_symbol=="NF1"& State=="MES-like")) %>% 
  mutate(alt = recode(snv, `1` = "mut",
                      `0` = "no mut"),
         genetic_alt = paste0(gene_symbol, " ", alt))

genetic_alt_order <- c("TERT no mut",
                       "TERT mut",
                       "TP53 no mut",
                       "TP53 mut", 
                       "NF1 no mut",
                       "NF1 mut")
care_all_snv_selected$genetic_alt <- factor(care_all_snv_selected$genetic_alt, levels=genetic_alt_order)
care_all_snv_selected$State <- factor(care_all_snv_selected$State, levels=c("GPC-like", "NEU-like", "MES-like", "Cycling"))
care_all_snv_selected$State <- droplevels(care_all_snv_selected$State)

res <- care_all_snv_selected %>% 
  filter(State%in%c("GPC-like", "MES-like", "NEU-like", "Cycling"))

# Inspect Wilcoxon rank sum test results - GPC-like (TERTp), NEU-like (TP53), MES-like (NF1), and Cycling (TP53).
ggplot(res, aes(x = as.factor(snv), y = freq*100)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox", label="p.format") +
  facet_grid(.~State, scales="free", space="free") +
  plot_theme +
  stat_n_text()

# CNV
cnv_state_long_hq_selected <- cnv_state_long_hq %>% 
  filter((gene_symbol=="MDM2" & State=="Cycling") | 
           (gene_symbol=="EGFR" & State=="GPC-like") | 
           (gene_symbol=="PDGFRA" & State=="OPC-like") |
           (gene_symbol=="CDK4"& State=="NPC-like") |
           (gene_symbol=="EGFR" & State=="AC-like") |
           (gene_symbol=="CDKN2A" & State=="GPC-like") |
           (gene_symbol=="7p" & State=="Hypoxia")) %>% 
  mutate(alt = recode(cnv_status, `1` = "amp",
                      `0` = "no amp"),
         genetic_alt = paste0(gene_symbol, " ", alt)) %>% 
  filter(!is.na(cnv_status))
cnv_state_long_hq_selected$genetic_alt <- gsub("CDKN2A no amp", "CDKN2A no del", cnv_state_long_hq_selected$genetic_alt) 
cnv_state_long_hq_selected$genetic_alt <- gsub("CDKN2A amp", "CDKN2A del", cnv_state_long_hq_selected$genetic_alt) 

genetic_alt_order <- c("EGFR no amp",
                       "EGFR amp",
                       "MDM2 no amp",
                       "MDM2 amp", 
                       "CDK4 no amp",
                       "CDK4 amp", 
                       "PDGFRA no amp",
                       "PDGFRA amp",
                       "CDKN2A no del",
                       "CDKN2A del",
                       "7p no amp",
                       "7p amp")
cnv_state_long_hq_selected$genetic_alt <- factor(cnv_state_long_hq_selected$genetic_alt, levels = genetic_alt_order)
cnv_state_long_hq_selected$State <- factor(cnv_state_long_hq_selected$State, levels = c("GPC-like", "AC-like", "Hypoxia", "Cycling", "NPC-like", "OPC-like"))

res_cnv <- cnv_state_long_hq_selected %>% 
  filter(State%in%c("GPC-like", "AC-like", "Hypoxia", "Cycling", "NPC-like", "OPC-like"))

# Examine uncorrected association between chromosome 7p amp and Hypoxia
ggplot(res_cnv %>% 
         filter(genetic_alt%in%c("7p no amp", "7p amp")), aes(x = as.factor(cnv_status), y = freq*100)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox", label="p.format") +
  facet_grid(.~State, scales="free", space="free") +
  plot_theme + 
  stat_n_text() # Note that some samples were excluded because of discontiguous segments

# Examine association between individual oncogene amplifications/deletions and malignant state abundance - unadjusted for multiple hypothesis testing
# NPC-like (CDK4), OPC-like (PDGFRA), and GPC-like (CDKN2A)
ggplot(res_cnv %>% 
         filter(genetic_alt%in%c("CDK4 no amp", "CDK4 amp", "PDGFRA no amp", "PDGFRA amp", "CDKN2A no del", "CDKN2A del")), aes(x = as.factor(cnv_status), y = freq*100)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox", label="p.format") +
  facet_grid(.~State, scales="free", space="free") +
  plot_theme + 
  stat_n_text()

# GPC-like (EGFR), AC-like (EGFR), and Cycling (MDM2) 
ggplot(res_cnv %>% 
         filter(genetic_alt%in%c("EGFR no amp", "EGFR amp", "MDM2 no amp", "MDM2 amp")), aes(x = as.factor(cnv_status), y = freq*100)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox", label="p.format") +
  facet_grid(.~State, scales="free", space="free") +
  plot_theme +
  stat_n_text()

# Combine both SNV and CNV into the same data.frame
res_snv_filtered <- res %>% 
  dplyr::select(ID, State, freq, genetic_alt)

res_cnv_filtered <- res_cnv %>% 
  dplyr::select(ID, State, freq, genetic_alt)

cnv_snv_comb <- bind_rows(res_snv_filtered, res_cnv_filtered)
cnv_snv_comb$State <- factor(cnv_snv_comb$State, levels=c("NPC-like", "OPC-like",  "MES-like", "AC-like", "GPC-like", "Hypoxia", "NEU-like", "Cycling"))
cnv_snv_comb$freq <- cnv_snv_comb$freq*100
jitter <- position_jitter(width = .15)

pdf("figures/nomura_fig5a_full_boxplot_points.pdf", width=9, height=5, useDingbats = FALSE, bg = "transparent")
ggplot(cnv_snv_comb, aes(x = genetic_alt, y = freq, fill=genetic_alt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = jitter, alpha = 0.25, size=0.75) +
  labs(x = "",
       y = "Relative malignant\ncell percentage (%)",
       fill = "Gene alteration") +
  facet_grid(.~State, scales="free", space="free") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("TERT no mut" = "gray80", 
                             "TERT mut" ="#33A02C",
                             "TP53 no mut" = "gray80", 
                             "TP53 mut" ="#33A02C",
                             "NF1 no mut" = "gray80", 
                             "NF1 mut" ="#33A02C", 
                             "EGFR amp" = "#E31A1C", 
                             "EGFR no amp" ="gray80",
                             "CDKN2A del" = "#1F78B4", 
                             "CDKN2A no del" ="gray80",
                             "MDM2 amp" = "#E31A1C", 
                             "MDM2 no amp" ="gray80",
                             "CDK4 amp" = "#E31A1C", 
                             "CDK4 no amp" ="gray80",
                             "PDGFRA amp" = "#E31A1C", 
                             "PDGFRA no amp" ="gray80",
                             "7p amp" = "#E31A1C", 
                             "7p no amp" ="gray80")) +
  theme(strip.text.x = element_text(margin = margin(t = 15, unit = "pt"))) +
  guides(fill=FALSE) +
  stat_n_text() 
dev.off()

# NOTE: All gene-level cna analyses here should be n = 108, all gene snv_level analyses should be n = 87, and chr7 should be 91 due to disrupted segments that made it difficult to assign CN status.

