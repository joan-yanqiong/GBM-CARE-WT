#################################
## Title: Figure 5a-e in Spitzer et al - longitudinal genetic alterations in cohort and associated gene expression changes
## Date: 2025.02.18
## Author: Kevin Johnson
## Description: Visualize longitudinal changes in mutational patterns and cell/expression changes 
#################################

# Necessary packages
library(tidyverse) # v1.3.1
library(ggpubr)
library(EnvStats)
library(viridis)
library(scales)
source("R/plot_theme.R")

# File containing information to link snrna and bulk dna identifiers.
linker_file_samples <- read.table("data/care_wt_synapse_dna_snrna_id_linker.txt", sep="\t", header = TRUE)
linker_file_cases <- linker_file_samples %>% 
  mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  dplyr::select(Patient, case_barcode) %>% 
  distinct()

### ### ### ### ### ###
### EDF6a - longitudinal changes in mutation burden
### ### ### ### ### ###
# Tumor mutation burden metrics both for matched blood and tumor-only calls (WARNING: tumor-only burden is inflated and only qualitatively used to assess relative changes).
longitudinal_mb_matched_tonly <- read.delim("data/longitudinal_mutation_burden.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Inspect normal blood and tumor matched longitudinal mutation burden changes.
matched_mb_plot <- ggplot(longitudinal_mb_matched_tonly %>% 
                          filter(status=="matched"), aes(x = timepoint, y = mf, fill=timepoint)) + 
  geom_hline(yintercept = 10, linetype=2, alpha=0.6) +
  geom_boxplot() +
  scale_y_log10() +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 3) +
  labs(x="", y="Mutation per megabase") +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  theme(strip.background = element_blank()) +
  stat_n_text() 

pdf(file = "figures/spitzer_edf6a_longitudinal_mb.pdf", height = 4, width = 4, useDingbats = FALSE)
matched_mb_plot
dev.off()

# There are two tumors with high mutation burden at initial time point with no prior surgery/treatment.
denovo_hm_cases <- longitudinal_mb_matched_tonly %>% 
  filter(status=="matched", timepoint=="T1", mf > 10 )
# Using the GLASS threshold for mutation burden c
hm_matched_cases <- longitudinal_mb_matched_tonly %>% 
  filter(status=="matched", timepoint=="T2", mf > 10 )
# Setting a threshold based on the mutation burden plot above for tumor-only. One sample identified with an 18x increase in mutation burden.
hm_tonly_cases <- longitudinal_mb_matched_tonly %>% 
  filter(status=="tumor-only", timepoint=="T2", mf > 50 )

# Define samples that have acquired hypermutation versus de novo hypermutation. De novo HM samples had no evidence of alkylating agent signature.
acquired_hypermutation <- hm_matched_cases %>% 
  bind_rows(hm_tonly_cases) %>% 
  filter(!case_barcode%in%denovo_hm_cases$case_barcode)

### ### ### ### ### ###
### EDF6b - De novo hypermutation mutational signature samples
### ### ### ### ### ###
# Single Base Substitution frequency for previously defined mutational signatures
cosmic_sbs_v3_care <- read.delim("data/cosmic_sbs_frequency.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

denovo_prop <- cosmic_sbs_v3_care %>% 
  mutate(case_barcode = substr(sample_barcode, 1, 12)) %>% 
  filter(case_barcode%in%c("CARE-TO-SM06", "CARE-TK-TK14")) %>% 
  mutate(signatures = recode(SBS.Sig.max, `SBS1` = "SBS1 (clock-like)",
                             `SBS5` = "SBS5 (clock-like)",
                             `SBS6` = "SBS6 (dMMR)",
                             `SBS8` = "SBS8 (Unknown)",
                             `SBS9` = "SBS9 (Pol eta)",
                             `SBS11` = "SBS11 (TMZ)",
                             `SBS12` = "SBS12 (Unknown)",
                             `SBS21` = "SBS21 (dMMR)",
                             `SBS30` = "SBS30 (dBER)",
                             `SBS39` = "SBS39 (Unknown)")) 
sig_levels <- c("SBS1 (clock-like)", "SBS5 (clock-like)", "SBS6 (dMMR)", "SBS8 (Unknown)", "SBS9 (Pol eta)",
                "SBS11 (TMZ)", "SBS12 (Unknown)", "SBS21 (dMMR)", "SBS30 (dBER)", "SBS39 (Unknown)")
denovo_prop$signatures <- factor(denovo_prop$signatures, levels=sig_levels)

pdf("figures/spitzer_edf6b_de_novo_hm.pdf", width = 5.5, height = 5, useDingbats = FALSE)
ggplot(denovo_prop, aes(x=ID, y=freq*100, fill=as.factor(signatures))) +
  geom_bar(stat="identity", position = "stack") +
  labs(x="De novo hypermutation samples", y="Relative percentage of mutations\n COSMICv3 SBS signatures (%)", fill="COSMIC V3 SBS") +
  plot_theme
dev.off()

### ### ### ### ### ###
### EDF6c - Glioma Longitudinal AnalySiS (GLASS) selected longitudinal mutational signatures
### ### ### ### ### ###
glass_cosmic_sbs_v3 <- read.delim("data/glass_cosmicv3_sbs_freq_n138.txt", sep="\t", header = TRUE)

pdf("figures/spitzer_edf6c_glass_mut_signatures.pdf", width = 6, height = 5, useDingbats = FALSE)
ggplot(glass_cosmic_sbs_v3 %>% 
         filter(SBS.Sig.max%in%c("SBS1", "SBS11", "SBS21")), aes(x = Timepoint, y = freq*100)) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=2) +
  scale_linetype_manual(values="dashed") +
  geom_boxplot(aes(fill=Timepoint)) +
  scale_fill_manual(values=c("T1"= "#66C2A5", "T2"="#FC8D62")) +
  plot_theme +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 3, label="p.format") +
  labs(x="GLASS - IDHwt (n = 138)", y="Relative proportion of mutations\n COSMICv3 SBS signatures") +
  theme(strip.background = element_blank()) +
  guides(fill=FALSE) +
  facet_grid(.~SBS.Sig.max, scales="free") +
  stat_n_text() 
dev.off()

### ### ### ### ### ###
### EDF6d - Longitudinal changes in small insertions
### ### ### ### ### ###
# Approximate breakdown of mutations and relative burden by variant type (deletions, insertions, etc)
longitudinal_variant <- read.delim("data/variant_type_longitudinal_breakdown.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Small insertion burden - for longitudinal analysis
longitudinal_ins <- longitudinal_variant %>% 
  filter(variant_type=="INS") %>% 
  # Remove the de novo hypermutator with many events at both timepoints
  filter(case_barcode!="CARE-TO-SM06") %>% 
  dplyr::select(case_barcode:tumor_barcode_b, variant_type, mf_initial= mf_is, mf_recurrence = mf_rs) %>% 
  pivot_longer(cols= c(mf_initial:mf_recurrence),
               names_to = "timepoint",
               values_to = "mf") %>% 
  mutate(timepoint = recode(timepoint, `mf_initial` = "T1",
                            `mf_recurrence` = "T2"),
         paired = rep(1:(n()/2),each=2))

small_ins_plot <- ggplot(longitudinal_ins %>% 
              # Exclude treatment-assoc. HM cases, which have many small insertions. If included, analysis remains statistically insignificant (p = 0.6).
              filter(!tumor_barcode_b%in%acquired_hypermutation$tumor_barcode_b), aes(x = timepoint, y = mf, fill=timepoint)) + 
  geom_line(aes(group=paired), color="gray70", linetype=2) +
  geom_boxplot() +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size =4, label="p.format") +
  labs(x="Time point", y="Small insertion\nper megabase") +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  theme(strip.background = element_blank()) +
  stat_n_text()

pdf(file = "figures/spitzer_edf6d_small_insertion.pdf", height = 4, width = 4, useDingbats = FALSE)
small_ins_plot
dev.off()


### ### ### ### ### ###
### EDF6e - Small deletion burden detected in recurrence versus delta (T2-T1) Hypoxia malignant state abundance
### ### ### ### ### ###
# Malignant cell annotation used throughout the manuscript.
malignant <- readRDS("data/malignant_meta_data_2025_01_08.RDS")

# Combine malignant cell proportions + cycling percentage.
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

malignant_freq = malignant %>%
  group_by(ID, State) %>%
  mutate(State = paste0(State, "-like")) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(ID, State,
           fill = list(counts = 0, freq = 0)) %>%
  mutate(Population = "malignant") %>% 
  # Append cycling percentage
  bind_rows(malignant_freq_cc) %>% 
  mutate(State = recode(State, `Hypoxia-like` = "Hypoxia",
                        `Neuron-like` = "NEU-like",
                        `Stress-like` = "Stress",
                        `Unresolved-like` = "Unresolved"))

malignant_cell_number_keep <- malignant %>% 
  group_by(ID) %>% 
  summarise(counts = n()) %>% 
  mutate(Patient =  sapply(strsplit(ID, "T"), "[[", 1),
         Timepoint= paste0("T", sapply(strsplit(ID, "T"), "[[", 2))) %>% 
  # Set the criteria for minimum number of malignant cells per sample
  filter(counts > 49, Timepoint!="T3") 

## Filter to those DNAseq files that pass cnv filters AND have sufficient malignant cell number.
malignant_freq_hq_wide <- malignant_freq %>% 
  filter(ID%in%malignant_cell_number_keep$ID) %>% 
  filter(Population%in%c("malignant", "cell_cycle")) %>% 
  mutate(Patient =  sapply(strsplit(ID, "T"), "[[", 1),
         Timepoint= paste0("T", sapply(strsplit(ID, "T"), "[[", 2))) %>% 
  dplyr::select(Patient, State, Timepoint, freq) %>%
  pivot_wider(names_from = Timepoint, values_from = freq) %>% 
  filter(!is.na(T1), !is.na(T2)) 

malignant_freq_hq_hypoxia <- malignant_freq_hq_wide %>% 
  filter(State == "Hypoxia") 

longitudinal_del <- longitudinal_variant %>% 
  filter(variant_type=="DEL") %>% 
  # Remove the de novo ultra-hypermutator with many events at both time points and defective mismatch repair (via SBS). The del frequency is about 6x higher in this sample than next highest sample - likely due to a germline event.
  filter(case_barcode!="CARE-TO-SM06") %>% 
  inner_join(linker_file_cases, by="case_barcode") %>% 
    filter(!case_barcode%in%acquired_hypermutation$case_barcode)

malignant_freq_hq_hypoxia_del <- malignant_freq_hq_hypoxia %>% 
  mutate(delta_hypoxia = T2-T1) %>% 
  inner_join(longitudinal_del, by="Patient")

n_distinct(malignant_freq_hq_hypoxia_del$Patient)
cor.test(malignant_freq_hq_hypoxia_del$delta_hypoxia, malignant_freq_hq_hypoxia_del$mf_r, method="p")

pdf("figures/spitzer_edf6e_smalldel_min50_n34.pdf", width = 5, height = 4, useDingbats = FALSE)
ggplot(malignant_freq_hq_hypoxia_del, aes(x= mf_r, y=delta_hypoxia*100)) +
  geom_point() +
  stat_cor(method="pearson", cor.coef.name = "R") +
  geom_smooth(method=lm, se = FALSE) +
  labs(x= "Small deletions acquired\nat recurrence (del Mb-1)", y="T2-T1 Hypoxia malignant\nstate abundance (%)") +
  plot_theme
dev.off()

### ### ### ### ### ###
### EDF6f - Longitudinal change in Hypoxia and MES-like for samples with increased SBS21 without acquired small deletions
### ### ### ### ### ###
sbs21_df <- read.delim("data/top_sbs21_samples.txt", header = TRUE, sep = "\t")
small_deletion_cases_only <- read.delim("data/small_deletion_cases.txt", header = TRUE, sep = "\t")

malignant_freq_hq <- malignant_freq %>% 
  filter(ID%in%malignant_cell_number_keep$ID) %>% 
  filter(Population%in%c("malignant", "cell_cycle")) %>% 
  mutate(Patient =  sapply(strsplit(ID, "T"), "[[", 1),
         Timepoint= paste0("T", sapply(strsplit(ID, "T"), "[[", 2))) %>% 
  filter(State%in%c("MES-like", "Hypoxia"))

malignant_freq_hq_annot <- malignant_freq_hq %>% 
  # Prioritize small deletion assignment
  mutate(phenotype = ifelse(Patient%in%unique(small_deletion_cases_only$Patient), "Small del", ifelse(Patient%in%unique(sbs21_df$Patient), "SBS21 increase", "Neither")),
         SBS21_phenotype = ifelse(Patient%in%unique(sbs21_df$Patient), "SBS21 increase", "Other"),
         Smalldel_phenotype = ifelse(Patient%in%unique(small_deletion_cases_only$Patient), "Small del", "Other")) %>% 
  filter(State%in%c("MES-like", "Hypoxia"),
         phenotype%in%c("SBS21 increase", "Small del")) %>% 
  filter(Patient%in%malignant_freq_hq_wide$Patient)

pdf(file = "figures/spitzer_edf6f_sbs21_state_changes.pdf", height = 4, width = 4, useDingbats = FALSE)
ggplot(malignant_freq_hq_annot %>% 
         filter(phenotype=="SBS21 increase", State%in%c("Hypoxia", "MES-like")), aes(x = Timepoint, y = freq*100)) + 
  geom_line(aes(group=Patient), color="gray70", linetype=2, alpha =0.5) +
  scale_color_manual(values=c("SBS21 increase"= "#FB6A4A", "Small del"="#A50F15")) +
  geom_point(aes(color=phenotype)) +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 3, label="p.format") +
  labs(x="Acquired SBS21 phenotype", y="Relative malignant cell abundance (%)") +
  theme(strip.background = element_blank()) +
  guides(color=FALSE) +
  facet_grid(.~State, scales="free") +
  plot_theme +
  stat_n_text()
dev.off()

# Significant increase when combining small deletion OR SBS21 increases: Hypoxia p = 0.00024 and MES-like p = 0.0012
ggplot(malignant_freq_hq_annot %>% 
         filter(State%in%c("Hypoxia", "MES-like")), aes(x = Timepoint, y = freq*100)) + 
  geom_line(aes(group=Patient), color="gray70", linetype=2, alpha =0.5) +
  scale_color_manual(values=c("SBS21 increase"= "#FB6A4A", "Small del"="#A50F15")) +
  geom_point(aes(color=phenotype)) +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 3, label="p.format") +
  labs(x="Acquired Small del. or SBS21 phenotype", y="Relative malignant cell abundance (%)") +
  theme(strip.background = element_blank()) +
  guides(color=FALSE) +
  facet_grid(.~State, scales="free") +
  plot_theme +
  stat_n_text()

### ### ### ### ### ###
### EDF6g - GLASS longitudinal changes in MES ssGSEA values for tumors with small deletion burden increases
### ### ### ### ### ###
# GLASS bulk RNA sequencing samples with both DNA and RNA that were classified as IDH-wildtype were scored for the MES-like (CARE) metaprogram
glass_dna_rna_smalldel_mes <- read.delim("data/glass_longitudinal_mes_ssgsea.txt", sep = "\t", header = TRUE)

pdf("figures/spitzer_edf6g_glass_mes.pdf", width = 3.5, height = 4, useDingbats = FALSE)
ggplot(glass_dna_rna_smalldel_mes %>% filter(type=="Acquired small del."), aes(x = timepoint, y = mes_ss_gsea,  fill=timepoint)) + 
  geom_line(aes(group=paired), color="gray70", linetype=2) +
  geom_boxplot() +
  scale_fill_manual(values=c("T1"= "#66C2A5", "T2"="#FC8D62")) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.format", method = "wilcox", paired = TRUE) +
  labs(y="Bulk RNA ssGSEA\nsnRNA MES metaprogram", x = "GLASS - longitudinal pairs") +
  theme(strip.background = element_blank()) +
  stat_n_text() +
  facet_grid(.~type, scales="free")
dev.off()

### ### ### ### ### ###
### EDF6h - Longitudinal transcriptomic distance (snRNA) versus genetic distance (mutations only)
### ### ### ### ### ###
# Require both time points to have at least 50 malignant cells. It ends up being n = 50 patients with at least 50 cells at both time points.
longitudinal_malignant_cell_analysis <- malignant %>% 
  group_by(ID) %>% 
  summarise(counts = n()) %>% 
  mutate(Patient =  sapply(strsplit(ID, "T"), "[[", 1),
         Timepoint= paste0("T", sapply(strsplit(ID, "T"), "[[", 2))) %>% 
  filter(counts > 49) %>% 
  dplyr::select(Patient, Timepoint, counts) %>% 
  pivot_wider(names_from = Timepoint, values_from = counts) %>% 
  filter(!is.na(T1), !is.na(T2))

# Load the transcriptomic distance metric between two samples - snRNA.
trans_distance_metric <- read.delim("data/transcriptomic_distance_tbl.csv", header = TRUE, sep = ",")

# Filter to the longitudinal set with minimum number of malignant cells at both time points.
trans_distance_metric_filt <- trans_distance_metric %>% 
  filter(Patient%in%unique(longitudinal_malignant_cell_analysis$Patient))

# Load the proportional change in CNA across CN segments. Intended to capture genetic distance from a CNA perspective.
gatk_seg_diff_prop <- read.delim("data/bulk_gatk_seg_diff_prop_n46.txt", header = TRUE, sep = "\t")

# Enumerating mutations where there is coverage of at least 15x to determine relative longitudinal mutation fraction.
mut_freq_prop_case <- readRDS("data/genetic_patient_mutation_fraction.RDS")
# Shared mutations
mut_freq_prop_shared <- mut_freq_prop_case %>% 
  filter(mutation_type=="S")


# A total of 38 cases that matches what is presented for transcriptional distance vs genetic distance.
genetic_transcript_dist_care <- mut_freq_prop_shared %>% 
  inner_join(trans_distance_metric_filt, by=c("patient"="Patient")) %>% 
  inner_join(gatk_seg_diff_prop, by=c("case_barcode")) %>% 
  mutate(Hypermutant = ifelse(case_barcode%in%hm_matched_cases$case_barcode, "HM", "nonHM"),
         mut_distance = 1-mutation_percent,
         genetic_distance = (mut_distance+prop_change)/2,
         scaled_genetic_distance = rescale(genetic_distance),
         scaled_snv_distance = rescale(mut_distance),
         scaled_cnv_distance = rescale(prop_change))

# There are 38 patients in this analysis.
n_distinct(genetic_transcript_dist_care$case_barcode)
# Pearson correlation coefficient= 0.382 and p-value = 0.0179.
cor.test(genetic_transcript_dist_care$scaled_snv_distance, genetic_transcript_dist_care$Dist01, method = "p")

pdf("figures/spitzer_edf6h_mutation_vs_transcriptomic_distance.pdf", width = 5, height = 4, useDingbats = FALSE)
ggplot(genetic_transcript_dist_care, aes(x=scaled_snv_distance, y=Dist01)) +
  stat_cor(method = "pearson", cor.coef.name = "R") +
  geom_point() +
  plot_theme +
  labs(x = "Genetic distance T1-T2\n(SNV diff.)", y="Transcriptional distance T1-T2\n(Euclidean distance metric)") +
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "bottom")
dev.off()



### ### ### ### ### ###
### EDF6i - Longitudinal transcriptomic distance (snRNA) versus genetic distance (CNAs only)
### ### ### ### ### ###
# Pearson correlation coefficient= 0.305334 and p-value = 0.0623.
cor.test(genetic_transcript_dist_care$scaled_cnv_distance, genetic_transcript_dist_care$Dist01, method = "p")

pdf("figures/spitzer_edf6i_cna_vs_transcriptomic_distance.pdf", width = 5, height = 4, useDingbats = FALSE)
ggplot(genetic_transcript_dist_care, aes(x=scaled_cnv_distance, y=Dist01)) +
  stat_cor(method = "pearson", cor.coef.name = "R") +
  geom_point() +
  plot_theme +
  labs(x = "Genetic distance T1-T2\n(CNA diff.)", y="Transcriptional distance T1-T2\n(Euclidean distance metric)") +
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "bottom")
dev.off()


### END ###