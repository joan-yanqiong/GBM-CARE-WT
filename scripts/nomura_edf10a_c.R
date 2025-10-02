#################################
## Title: Extended Data Figure 10 panels a-c in Nomura et al - genetic alterations and malignant cell state abundance
## Date: 2025.02.20
## Author: Kevin Johnson
## Description: Common GBM genetic alteration and malignant cell state abundance analyses
#################################

# Necessary packages
library(tidyverse) # v1.3.1
library(ggpubr) # v0.6.0
library(EnvStats) # v2.8
library(broom)
library(openxlsx)
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


## Examine malignant frequency. Tabulate cell cycle first.
malignant$Patient <-  sapply(strsplit(malignant$ID, "T"), "[[", 1)
malignant$Timepoint <-  paste0("T", sapply(strsplit(malignant$ID, "T"), "[[", 2))
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

# Set threshold for samples that should be kept in an analysis at more than 50 malignant cells.
malignant_cell_number_keep <- malignant_cell_number %>% 
  filter(counts >50)

# Creating a data frame with "high-quality" malignant cells. This is to reduce noise associated with samples that only have a few malignant cells.
# This goes from 120 samples (1 dropped due to no malignant cells) to 111 samples with at greater than 50 malignant cells.
all_malignant_freq_hq <- all_cell_freq %>% 
  filter(ID%in%malignant_cell_number_keep$ID) %>% 
  mutate(Patient = sapply(strsplit(ID, "T"), "[[", 1)) %>% 
  mutate(State = recode(State, `Hypoxia-like` = "Hypoxia",
                        `Neuron-like` = "NEU-like",
                        `Stress-like` = "Stress",
                        `Unresolved-like` = "Unresolved"))

### ### ### ### ### ###
### CNA arm heatmap
### ### ### ### ### ###
# First assess whether there are any significant differences via a wilcoxon test.
wilcox_arm_res =  gatk_cnv_by_arm_care %>% 
  dplyr::select(ID, aliquot_barcode, arm, arm_call) %>% 
  # Selecting chromosome arms frequently altered in glioma. Need to select a specific direction for alterations (i.e., gain or loss).
  filter(arm%in%c("7p","7q", "19p", "19q", "20p", "20q", "10p", "10q", "13q", "14q", "22q")) %>% 
  mutate(cnv_status = ifelse(arm %in% c("7p","7q", "19p", "19q", "20p", "20q") & arm_call == 1, 1,
                             ifelse(arm %in% c("10p", "10q", "13q", "14q", "22q") & arm_call == -1, 1, 0))) %>% 
  inner_join(all_malignant_freq_hq, by=c("ID")) %>% 
  dplyr::select(c(ID, State, freq, arm, cnv_status)) %>% 
  group_nest(State, arm) %>%
  mutate(wilcox_test = map(data, ~tidy(wilcox.test(freq ~ cnv_status, data = .x)))) %>%
  unnest(wilcox_test)

# Format data to understand the direction of association.
arm_state_long <- gatk_cnv_by_arm_care %>% 
  dplyr::select(ID, aliquot_barcode, arm, arm_call) %>% 
  filter(arm%in%c("7p","7q", "19p", "19q", "20p", "20q", "10p", "10q", "13q", "14q", "22q")) %>% 
  mutate(cnv_status = ifelse(arm %in% c("7p","7q", "19p", "19q", "20p", "20q") & arm_call == 1, 1,
                             ifelse(arm %in% c("10p", "10q", "13q", "14q", "22q") & arm_call == -1, 1, 0))) %>% 
  inner_join(all_malignant_freq_hq, by=c("ID")) 

arm_state_wide <- arm_state_long %>% 
  group_by(State, arm, cnv_status) %>% 
  summarise(median_state_freq = median(freq)) %>% 
  filter(!is.na(cnv_status)) %>% 
  pivot_wider(names_from = cnv_status, values_from = median_state_freq) %>% 
  dplyr::select(State, arm, non_alt = `0`, alt = `1`) %>% 
  mutate(delta = alt-non_alt,
         direction = ifelse(delta>0, "increased_w_alt", "decreased_w_alt")) %>% 
  ungroup()

# Add information about median value and perform p-value adjustment for multiple hypothesis testing.
wilcox_arm_res_sig <- wilcox_arm_res %>% 
  inner_join(arm_state_wide, by=c("State", "arm")) %>% 
  mutate(adj.p.value = p.adjust(p.value, method = "fdr")) %>% 
  # Mainly, we are interest in increases in cell state abundance due to an alteration.
  filter(direction == "increased_w_alt") %>% 
  dplyr::select(State, arm, direction, p.value, adj.p.value) %>% 
  # Restricting to the high frequency events.
  mutate(arm_alt = ifelse(arm %in% c("7p","7q", "19p", "19q", "20p", "20q"), paste0(arm, " amp"),
                          ifelse(arm %in% c("10p", "10q", "13q", "14q", "22q"),  paste0(arm, " del"), arm)))

malignant_state_order <- c("Cilia-like" , 
                           "AC-like",
                           "Hypoxia",
                           "MES-like",
                           "Stress",
                           "GPC-like",
                           "OPC-like",
                           "NPC-like",
                           "NEU-like",
                           "Cycling",
                           "Unresolved")
wilcox_arm_res_sig$State <- factor(wilcox_arm_res_sig$State, levels=malignant_state_order)
wilcox_arm_res_sig$arm_alt <- factor(wilcox_arm_res_sig$arm_alt, levels=rev(c("7p amp","7q amp", "19p amp", "19q amp", "20p amp", "20q amp", "10p del", "10q del", "13q del", "14q del", "22q del")))

# Note: Less confidence around 14q del and Cilia-like because Cilia-like abundance is driven by a handful of samples

pdf("figures/arm_cnv_malignant_state_all.pdf", width=8, height=5, useDingbats = FALSE, bg = "transparent")
ggplot(wilcox_arm_res_sig %>% 
         filter(arm%in%c("7p", "7q", "19p", "14q", "13q", "10q", "10p")), aes(x = State, y = arm_alt, fill=-log10(adj.p.value))) +
  geom_tile() +
  scale_fill_gradient2("-log10(adj. P)\nsig. results", limits = c(1.3, max(-log10(wilcox_arm_res_sig$adj.p.value))), na.value = "white") +
  labs(x = "Malignant cell state", y = "Arm-level alteration") +
  theme_bw() +
  theme(axis.text.y = element_text(size=10), axis.text.x= element_text(size=10,angle=45,hjust=1),
        axis.title = element_text(size=10),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=10),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.15, 'inch')) +
  ggtitle("Genetic alterations associated with increases in malignant proportions")
dev.off()

### ### ### ### ### ###
### CNA heatmap
### ### ### ### ### ###
## Restrict to the genes previously associated with cell states in bulk TCGA scoring analyses and general GBM CNA driver events
wilcox_cnv_res = gene_copy_number_care %>% 
  dplyr::select(ID, aliquot_barcode, gene_symbol, hlvl_call) %>% 
  filter(gene_symbol%in%c("EGFR", "CDK4", "PDGFRA", "MDM2", "CDKN2A", "PTEN", "RB1")) %>% 
  mutate(cnv_status = ifelse(gene_symbol%in% c("EGFR", "CDK4", "PDGFRA", "MDM2") & hlvl_call == 2, 1,
                             # Selecting a more relaxed threshold of -1 and -2 for deletions due to the purity differences in tumors (e.g., a deep deletion may be called single copy loss in low purity tumor).
                             ifelse(gene_symbol%in%c("CDKN2A", "PTEN", "RB1") & hlvl_call%in%c(-1, -2), 1, 0))) %>% 
  inner_join(all_malignant_freq_hq, by=c("ID")) %>% 
  dplyr::select(c(ID, State, freq, gene_symbol, cnv_status)) %>% 
  group_nest(State, gene_symbol) %>%
  mutate(wilcox_test = map(data, ~tidy(wilcox.test(freq ~ cnv_status, data = .x)))) %>%
  unnest(wilcox_test)

## Restrict to the genes of interest:
gene_copy_number_care_long = gene_copy_number_care %>% 
  dplyr::select(ID, aliquot_barcode, gene_symbol, hlvl_call) %>% 
  filter(gene_symbol%in%c("EGFR", "CDK4", "PDGFRA", "MDM2", "CDKN2A", "PTEN", "RB1")) %>% 
  mutate(cnv_status = ifelse(gene_symbol%in% c("EGFR", "CDK4", "PDGFRA", "MDM2") & hlvl_call == 2, 1,
                             # Selecting a more relaxed threshold of -1 and -2 for deletions due to the purity differences in tumors (e.g., a deep deletion may be called single copy loss in low purity tumor).
                             ifelse(gene_symbol%in%c("CDKN2A", "PTEN", "RB1") & hlvl_call%in%c(-1, -2), 1, 0))) %>% 
  inner_join(all_malignant_freq_hq, by=c("ID"))

# Determine median levels:
cnv_state_wide <- gene_copy_number_care_long %>% 
  group_by(State, gene_symbol, cnv_status) %>% 
  summarise(median_state_freq = median(freq)) %>% 
  filter(!is.na(cnv_status)) %>% 
  pivot_wider(names_from = cnv_status, values_from = median_state_freq) %>% 
  dplyr::select(State, gene_symbol, non_alt = `0`, alt = `1`) %>% 
  mutate(delta = alt-non_alt,
         direction = ifelse(delta>0, "increased_w_alt", "decreased_w_alt")) %>% 
  ungroup()

# Inspect only those alterations that increase a specific malignant cell state
wilcox_cnv_res_sig <- wilcox_cnv_res %>% 
  inner_join(cnv_state_wide, by=c("State", "gene_symbol")) %>% 
  mutate(adj.p.value = p.adjust(p.value, method = "fdr")) %>% 
  filter(direction == "increased_w_alt") %>% 
  dplyr::select(State, gene_symbol, direction, p.value, adj.p.value) %>% 
  mutate(gene_alt = ifelse(gene_symbol%in% c("EGFR", "CDK4", "PDGFRA", "MDM2"), paste0(gene_symbol, " amp"),
                           ifelse(gene_symbol%in%c("CDKN2A", "PTEN", "RB1"),  paste0(gene_symbol, " del"), gene_symbol)))


# Re-order the factors
wilcox_cnv_res_sig$State <- factor(wilcox_cnv_res_sig$State, levels=malignant_state_order)
wilcox_cnv_res_sig$gene_alt <- factor(wilcox_cnv_res_sig$gene_alt, levels=rev(c("EGFR amp", "CDK4 amp", "PDGFRA amp", "MDM2 amp", "RB1 del","PTEN del", "CDKN2A del")))

# Note that RB1 and PTEN deletions are likely to be more reflective of chr10 and chr13 loss.
pdf("figures/gene_cnv_malignant_state_all.pdf", width=8, height=5, useDingbats = FALSE, bg = "transparent")
ggplot(wilcox_cnv_res_sig %>% 
         filter(gene_symbol%in%c("EGFR", "CDK4", "PDGFRA", "MDM2", "CDKN2A", "PTEN", "RB1")), aes(x = State, y = gene_alt, fill=-log10(adj.p.value))) +
  geom_tile() +
  scale_fill_gradient2("-log10(adj. P)\nsig. results", limits = c(1.3, max(-log10(wilcox_cnv_res_sig$adj.p.value))), na.value = "white") +
  labs(x = "Malignant cell state", y = "Gene-level CNV alteration") +
  theme_bw() +
  theme(axis.text.y = element_text(size=10), axis.text.x= element_text(size=10,angle=45,hjust=1),
        axis.title = element_text(size=10),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=10),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.15, 'inch')) +
  ggtitle("Genetic alterations associated with increases in malignant proportions")
dev.off()

# Unadjusted p-values. Among these nominally significant hits, those with previous bulk associations with CNAs (i.e., Neftel et al) are reported. Plus, new significant hits.
ggplot(wilcox_cnv_res_sig %>% 
         filter(gene_symbol%in%c("EGFR", "CDK4", "PDGFRA", "MDM2", "CDKN2A", "PTEN", "RB1")), aes(x = State, y = gene_alt, fill=-log10(p.value))) +
  geom_tile() +
  scale_fill_gradient2("-log10(adj. P)\nsig. results", limits = c(1.3, max(-log10(wilcox_cnv_res_sig$p.value))), na.value = "white") +
  labs(x = "Malignant cell state", y = "Gene-level CNV alteration") +
  theme_bw() +
  theme(axis.text.y = element_text(size=10), axis.text.x= element_text(size=10,angle=45,hjust=1),
        axis.title = element_text(size=10),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=10),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.15, 'inch')) +
  ggtitle("Genetic alterations associated with increases in malignant proportions")

### ### ### ### ### ###
### SNV heatmap
### ### ### ### ### ###
care_snv_long <- care_snv_df %>% 
  dplyr::select(ID, TERT:PIK3CA) %>% 
  pivot_longer(cols= c(TERT:PIK3CA),
               names_to = "gene_symbol",
               values_to = "snv") 

all_malignant_freq_hq_genetics_long <- all_malignant_freq_hq %>% 
  inner_join(care_snv_long, by="ID")

wilcox_test_res_snv = all_malignant_freq_hq_genetics_long %>% 
  dplyr::select(c(ID, State, freq, gene_symbol, snv)) %>% 
  group_nest(State, gene_symbol) %>%
  mutate(wilcox_test = map(data, ~tidy(wilcox.test(freq ~ snv, data = .x)))) %>%
  unnest(wilcox_test) 


# Determine median levels:
snv_state_wide <- all_malignant_freq_hq_genetics_long %>% 
  group_by(State, gene_symbol, snv) %>% 
  summarise(median_state_freq = median(freq)) %>% 
  pivot_wider(names_from = snv, values_from = median_state_freq) %>% 
  dplyr::select(State, gene_symbol, non_alt = `0`, alt = `1`) %>% 
  mutate(delta = alt-non_alt,
         direction = ifelse(delta>0, "increased_w_alt", "decreased_w_alt")) %>% 
  ungroup()

# Focus on those alterations that are associated with an increase in a specific cell state abundance.
wilcox_snv_res_sig <- wilcox_test_res_snv %>% 
  inner_join(snv_state_wide, by=c("State", "gene_symbol")) %>% 
  mutate(adj.p.value = p.adjust(p.value, method = "fdr")) %>% 
  filter(direction == "increased_w_alt") %>% 
  dplyr::select(State, gene_symbol, direction, p.value, adj.p.value) 

# Re-order the features
wilcox_snv_res_sig$State <- factor(wilcox_snv_res_sig$State, levels=malignant_state_order)

pdf("figures/gene_snv_malignant_state_all.pdf", width=8, height=5, useDingbats = FALSE, bg = "transparent")
ggplot(wilcox_snv_res_sig, aes(x = State, y = gene_symbol, fill=-log10(adj.p.value))) +
  geom_tile() +
  scale_fill_gradient2("-log10(adj. P)\nsig. results", limits = c(1.3, max(-log10(wilcox_cnv_res_sig$adj.p.value))), na.value = "white") +
  labs(x = "Malignant cell state", y = "Gene-level SNV alteration") +
  theme_bw() +
  theme(axis.text.y = element_text(size=10), axis.text.x= element_text(size=10,angle=45,hjust=1),
        axis.title = element_text(size=10),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=10),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.15, 'inch')) +
  ggtitle("Genetic alterations associated with increases in malignant proportions")
dev.off()

# Unadjusted significant p-values, these are the significant hits that get reported in Spitzer et al Figure 5a (e.g., NF1 and MES-like).
ggplot(wilcox_snv_res_sig, aes(x = State, y = gene_symbol, fill=-log10(p.value))) +
  geom_tile() +
  scale_fill_gradient2("-log10(adj. P)\nsig. results", limits = c(1.3, max(-log10(wilcox_cnv_res_sig$p.value))), na.value = "white") +
  labs(x = "Malignant cell state", y = "Gene-level SNV alteration") +
  theme_bw() +
  theme(axis.text.y = element_text(size=10), axis.text.x= element_text(size=10,angle=45,hjust=1),
        axis.title = element_text(size=10),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=10),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.15, 'inch')) +
  ggtitle("Genetic alterations associated with increases in malignant proportions")

wilcox_arm_res_sig_c <- wilcox_arm_res_sig %>% mutate(alt_type = "Chr. arm CNA", type = "Cell Abundance") %>%  dplyr::select(State, type, genetic_event = arm_alt, alt_type, p.value, adj.p.value)
wilcox_cnv_res_sig_c <- wilcox_cnv_res_sig %>% mutate(alt_type = "Gene CNA", type ="Cell Abundance")  %>%  dplyr::select(State, type, genetic_event = gene_alt, alt_type, p.value, adj.p.value)
wilcox_snv_res_sig_c  <- wilcox_snv_res_sig %>% mutate(alt_type = "SNV", type = "Cell Abundance") %>%  dplyr::select(State, type, genetic_event = gene_symbol, alt_type, p.value, adj.p.value)
wilcox_all_res_sig <- bind_rows(wilcox_arm_res_sig_c, wilcox_cnv_res_sig_c, wilcox_snv_res_sig_c)

wilcox_all_res_sig$State <- factor(wilcox_all_res_sig$State, levels=malignant_state_order)

# In an effort to reduce whitespace, only specific alterations are featured.
pdf("figures/nomura_edf10a_states_genetics_heatmap.pdf", width=6, height=5, useDingbats = FALSE, bg = "transparent")
ggplot(wilcox_all_res_sig %>% 
         filter(State%in%c("Hypoxia", "GPC-like", "NEU-like", "Cycling"), genetic_event%in%c("7p amp", "13q del", "10p del", 
                                                                                                     "EGFR amp", "MDM2 amp","PTEN del", "CDKN2A del",
                                                                                                     "TERT", "TP53")), aes(x = State, y = genetic_event, fill=-log10(adj.p.value))) +
  geom_tile() +
  scale_fill_gradient2("-log10(adj. P)\nsig. results", limits = c(1.3, max(-log10(wilcox_all_res_sig$adj.p.value))), na.value = "white") +
  labs(x = "Malignant cell state", y = "Genetic alteration") +
  theme_bw() +
  theme(axis.text.y = element_text(size=10), axis.text.x= element_text(size=10,angle=45,hjust=1),
        axis.title = element_text(size=10),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=10),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.15, 'inch')) +
  ggtitle("Genetic alt linked with increased state abundance") +
  facet_grid(alt_type~., scales = "free", space = "free")
dev.off()

### ### ### ### ###
### EDF10b - Any classic GBM genetic alteration (TERTp, EGFRamp, Chr7, Chr10) and GPC-like abundance
### ### ### ### ###
care_cnv_df <- gene_copy_number_care %>% 
  dplyr::select(ID, aliquot_barcode, gene_symbol, hlvl_call) %>% 
  filter(gene_symbol%in%c("EGFR", "CDK4", "PDGFRA", "MDM2", "CDKN2A", "PTEN", "RB1")) %>% 
  mutate(cnv_status = ifelse(gene_symbol%in% c("EGFR", "CDK4", "PDGFRA", "MDM2") & hlvl_call == 2, 1,
                             ifelse(gene_symbol%in%c("CDKN2A", "PTEN", "RB1") & hlvl_call%in%c(-1, -2), 1, 0))) %>% 
  dplyr::select(-hlvl_call) %>%
  mutate(gene_symbol = paste0(gene_symbol, "_cna")) %>% 
  pivot_wider(names_from = gene_symbol, values_from = cnv_status) 

# Arm-level 
care_arm_df <- gatk_cnv_by_arm_care %>% 
  dplyr::select(ID, aliquot_barcode, arm, arm_call) %>% 
  filter(arm%in%c("7p","7q", "19p", "19q", "20p", "20q", "10p", "10q", "13q", "14q", "22q")) %>% 
  mutate(cnv_status = ifelse(arm %in% c("7p","7q", "19p", "19q", "20p", "20q") & arm_call == 1, 1,
                             ifelse(arm %in% c("10p", "10q", "13q", "14q", "22q") & arm_call == -1, 1, 0))) %>% 
  dplyr::select(-arm_call) %>%
  pivot_wider(names_from = arm, values_from = cnv_status) 

care_all_genetics <- care_snv_df %>% 
  inner_join(care_cnv_df, by=c("ID","aliquot_barcode")) %>% 
  inner_join(care_arm_df, by=c("ID","aliquot_barcode")) 

care_all_genetics_classic <- care_all_genetics %>% 
  dplyr::select(ID:case_barcode, TERT, EGFR_cna, `CDKN2A_cna`,`7p`, `10p`) %>% 
  pivot_longer(cols= c(TERT:`10p`),
               names_to = "alt",
               values_to = "event_status")   %>% 
  group_by(ID, case_barcode) %>% 
  summarise(num_events = sum(event_status, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(genetic_class = ifelse(num_events<1, "no alt", ifelse(num_events>1 & num_events<3, "alt", "alt"))) 

# Most samples have multiple events
sum(care_all_genetics_classic$num_events>1)/95
sum(care_all_genetics_classic$num_events>=3)/95

care_all_genetics_classic_filt <- care_all_genetics_classic %>% 
  filter(num_events < 2) %>% 
  group_by(case_barcode) %>% 
  summarise(counts = n())

# There are 87 samples with both sufficient genetics and high quality malignant cells - due to chr7/chr10 missingness.
care_all_genetics_classic_malignant <- care_all_genetics_classic %>% 
  inner_join(all_malignant_freq_hq, by="ID")

care_all_genetics_classic_malignant$genetic_class <- factor(care_all_genetics_classic_malignant$genetic_class, levels=c("no alt", "alt"))
jitter <- position_jitter(width = .15)

pdf("figures/nomura_fig10b_multiple_classic_alt.pdf", width=1.5, height=2.2, useDingbats = FALSE, bg = "transparent")
ggplot(care_all_genetics_classic_malignant %>% 
         filter(State=="GPC-like"), aes(x = genetic_class, y = freq*100, fill=genetic_class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = jitter, alpha = 0.25) +
  labs(x = "",
       y = "GPC-like relative malignant\ncell abundance (%)",
       fill = "GBM State") +
  guides(fill=FALSE) +
  scale_fill_manual(values=c("alt" = "red", 
                             "no alt" ="gray80")) +
  theme_bw(base_size = 10) + theme(axis.title = element_text(size = 10),
                                   axis.text = element_text(size = 10),
                                   panel.background = element_rect(fill = "transparent"),
                                   axis.line = element_blank(),
                                   strip.background = element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), 
                                   panel.border = element_blank(),
                                   axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                   axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) +
  scale_x_discrete(labels=c('No alt', 'Alt')) +
  stat_compare_means(method="wilcox", label="p.format") +
  stat_n_text()
dev.off()

### ### ### ### ###
### EDF10c - Any classic GBM genetic alteration (TERTp, EGFRamp,Chr7,Chr10) and BP score
### ### ### ### ###
bp_class <- readWorkbook("data/nomura_supptables.xlsx", sheet = 5, startRow = 4, colNames = TRUE)

bp_class_long <- bp_class %>% 
  dplyr::select(ID=SampleID, BP_ECM:BP_Glial) %>% 
  pivot_longer(cols= c(BP_ECM:BP_Glial),
               names_to = "bp",
               values_to = "scores") %>%
  inner_join(care_all_genetics_classic, by="ID")

pdf("figures/nomura_edf10c_baseline_profile_scores.pdf", width=6, height=5, useDingbats = FALSE, bg = "transparent")
ggplot(bp_class_long, aes(x=as.factor(genetic_class), y=scores, fill=as.factor(genetic_class))) + 
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  stat_compare_means(method="wilcox.test", label="p.format") +
  stat_n_text() + 
  theme_bw() +
  scale_fill_manual(values=c("alt" = "red", 
                             "no alt" ="gray80")) +
  scale_x_discrete(labels=c('Alt', 'No alt')) +
  labs(x="Any EGFR/Chr7/Chr10/TERTp alt", y="Baseline Profile Score") +
  plot_theme +
  theme(legend.position="none") + 
  facet_grid(.~bp, scales="free", space="free") 
dev.off()

### END ###