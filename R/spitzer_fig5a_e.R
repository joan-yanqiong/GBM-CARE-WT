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
### Figure 5a - longitudinal differences in mutational processes/signatures
### ### ### ### ### ###
# Single Base Substitution frequency for previously defined mutational signatures
cosmic_sbs_v3_care <- read.delim("data/cosmic_sbs_frequency.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# There are a few T3 samples and a singleton T1 sample. Isolating the longitudinal T1 and T2 pairs.
cosmic_sbs_v3_care_pairs <- cosmic_sbs_v3_care %>% 
  dplyr::select(Patient, Timepoint) %>% 
  distinct() %>% 
  group_by(Patient) %>% 
  summarise(counts = n()) %>% 
  filter(counts > 1)

sbs_levels <- c("SBS1", "SBS5", "SBS6", "SBS7b", "SBS8", "SBS9", "SBS11", "SBS12", "SBS21", "SBS30", "SBS31", "SBS39")
cosmic_sbs_v3_care$SBS.Sig.max <- factor(cosmic_sbs_v3_care$SBS.Sig.max, levels=sbs_levels)

ggplot(cosmic_sbs_v3_care %>% filter(Patient%in%cosmic_sbs_v3_care_pairs$Patient, Timepoint!="T3"), aes(x = Timepoint, y = freq*100)) + 
  geom_line(aes(group=Patient), color="gray70", linetype=2) +
  geom_boxplot(aes(fill=SBS.Sig.max)) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 3, label="p.format") +
  labs(x="CARE - IDHwt pairs (n = 46)", y="Relative proportion of mutations\n COSMICv3 SBS signatures") +
  theme(strip.background = element_blank()) +
  guides(fill=FALSE) +
  facet_grid(.~SBS.Sig.max, scales="free") +
  stat_n_text() 

pdf("figures/spitzer_fig5a.pdf", width = 6, height = 5, useDingbats = FALSE)
ggplot(cosmic_sbs_v3_care %>% filter(Patient%in%cosmic_sbs_v3_care_pairs$Patient,
                                     SBS.Sig.max%in%c("SBS1", "SBS11", "SBS21"),
                                     Timepoint!="T3"), aes(x = Timepoint, y = freq*100)) + 
  geom_line(aes(group=Patient), color="gray70", linetype=2) +
  geom_boxplot(aes(fill=Timepoint)) +
  scale_fill_manual(values=c("T1"= "#66C2A5", "T2"="#FC8D62")) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 3, label="p.format") +
  labs(x="Tumor pairs (n = 46 pairs)", y="Relative proportion of mutations\n COSMICv3 SBS signatures") +
  theme(strip.background = element_blank()) +
  guides(fill=FALSE) +
  facet_grid(.~SBS.Sig.max, scales="free") +
  stat_n_text()
dev.off()

# Define which samples are in the top quartile for SBS21 change (i.e., SBS21 increases).
sbs21_wide_all <- cosmic_sbs_v3_care %>% 
  filter(Patient%in%cosmic_sbs_v3_care_pairs$Patient, SBS.Sig.max=="SBS21", Timepoint!="T3") %>% 
  dplyr::select(SBS.Sig.max, freq, Patient, Timepoint) %>% 
  pivot_wider(names_from = Timepoint, values_from = freq) %>% 
  mutate(sbs21_diff = T2-T1) 

summary(sbs21_wide_all$sbs21_diff)
summary(sbs21_wide_all$sbs21_diff)[[5]]
third_quartile <- summary(sbs21_wide_all$sbs21_diff)[[5]]

sbs21_wide <- sbs21_wide_all %>% 
  # Select samples based on the top quartile
  filter(sbs21_diff>third_quartile) %>%
  inner_join(linker_file_cases, by="Patient")

sbs21_wide_phenotype_out <- sbs21_wide %>% 
  mutate(Timepoint = "T2") %>% 
  mutate(ID = paste0(Patient, Timepoint)) %>% 
  dplyr::select(ID, Patient, Timepoint, Phenotype = SBS.Sig.max)

# save the top quartile samples with increases in SBS21.
write.table(sbs21_wide_phenotype_out, "data/top_sbs21_samples.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


### ### ### ### ### ###
### Figure 5b - longitudinal small deletion changes
### ### ### ### ### ###
# Tumor mutation burden metrics both for matched blood and tumor-only calls (WARNING: tumor-only burden is inflated and only qualitatively used to assess relative changes).
longitudinal_mb_matched_tonly <- read.delim("data/longitudinal_mutation_burden.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Inspect tumor-only longitudinal mutation burden changes.
tonly_mb_plot <- ggplot(longitudinal_mb_matched_tonly %>% 
               filter(status=="tumor-only"), aes(x = timepoint, y = mf, fill=timepoint)) + 
  geom_line(aes(group=paired), color="gray70", linetype=2) +
  geom_boxplot() +
  scale_y_log10() +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 3) +
  labs(x="", y="Mutation per megabase\n(tumor-only mutation calls)") +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  theme(strip.background = element_blank()) +
  stat_n_text() 
tonly_mb_plot

# There are two tumors with high mutation burden at initial time point with no prior surgery/treatment.
denovo_hm_cases <- longitudinal_mb_matched_tonly %>% 
  filter(status=="matched", timepoint=="T1", mf > 10 )
# Using the GLASS threshold for mutation burden c
hm_matched_cases <- longitudinal_mb_matched_tonly %>% 
  filter(status=="matched", timepoint=="T2", mf > 10 )
# All non-de novo hypermutant tumors with high mutation also have high SBS11, consistent with treatment-associated hypermutation. These were also the same when inspecting v2 signatures.
high_sbs11 <- cosmic_sbs_v3_care %>% 
  filter(SBS.Sig.max=="SBS11", freq > 0.4)
  # Setting a threshold based on the mutation burden plot above for tumor-only. One sample identified with an 18x increase in mutation burden.
hm_tonly_cases <- longitudinal_mb_matched_tonly %>% 
  filter(status=="tumor-only", timepoint=="T2", mf > 50 )

# Define samples that have acquired hypermutation versus de novo hypermutation. De novo HM samples had no evidence of alkylating agent signature.
acquired_hypermutation <- hm_matched_cases %>% 
  bind_rows(hm_tonly_cases) %>% 
  filter(!case_barcode%in%denovo_hm_cases$case_barcode)

# Breakdown by variant count.
longitudinal_variant <- read.delim("data/variant_type_longitudinal_breakdown.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Define acquired small deletion burden.
longitudinal_variant_del <- longitudinal_variant %>% 
  # Focus on small deletions.
  filter(variant_type=="DEL") %>% 
  mutate(small_del_group = ifelse(mf_r >=.2, 2, 1),
         mf_diff = mf_rs-mf_is) %>% 
  # Note that there are two requirements in this setting:
  # 1. A sample needs to have a small deletion mutation burden for its recurrence-only variants > 0.2 mutations/Mb.
  # 2. There needs to be a difference of 0.1 Mut/Mb between the recurrence and the initial tumor so that it is truly acquired.
  mutate(rt_scars = ifelse(small_del_group == 2 & mf_diff>0.1, "rt_scars @ T2", "no_rt_scars @ T2")) %>% 
  # Remove the de novo ultra-hypermutator with many events at both time points and defective mismatch repair (via SBS). The del frequency is about 6x higher in this sample than next highest sample - likely due to a germline event.
  filter(case_barcode!="CARE-TO-SM06") 

# Restrict to subjects of interest - 12 total samples, 10 of which are not acquired hypermutant tumors
small_del_samples <- longitudinal_variant_del %>% 
  filter(rt_scars=="rt_scars @ T2") %>% 
  dplyr::select(case_barcode:tumor_barcode_b) %>% 
  inner_join(linker_file_cases, by="case_barcode")

small_del_aliquots <- small_del_samples$tumor_barcode_b[which(!small_del_samples$tumor_barcode_b%in%acquired_hypermutation$tumor_barcode_b)]
small_del_cases <- substr(small_del_aliquots, 1, 12)

small_del_cases_df <- small_del_samples %>% 
  filter(!tumor_barcode_b%in%acquired_hypermutation$tumor_barcode_b)

# Write out result for other scripts to use small deletion burden classification.
write.table(small_del_cases_df, "data/small_deletion_cases.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


small_del_plot_input <- longitudinal_variant_del %>% 
  dplyr::select(case_barcode:tumor_barcode_b, variant_type, mf_initial= mf_is, mf_recurrence = mf_rs) %>% 
  pivot_longer(cols= c(mf_initial:mf_recurrence),
               names_to = "timepoint",
               values_to = "mf") %>% 
  mutate(timepoint = recode(timepoint, `mf_initial` = "T1",
                            `mf_recurrence` = "T2"),
         paired = rep(1:(n()/2),each=2)) %>% 
  mutate(smalldel = ifelse(case_barcode%in%small_del_cases, "smalldel", "other"))


# The total sample size here will be n = 42 because of removal of 3 acquired HM and the de novo hypermutant tumor with dMMR.
small_del_plot <- ggplot(small_del_plot_input %>% 
              filter(!tumor_barcode_b%in%acquired_hypermutation$tumor_barcode_b), aes(x = timepoint, y = mf, fill=timepoint)) + 
  stat_compare_means(method = "wilcox", paired = TRUE, size =4) +
  labs(x="", y="Small deletion\nper megabase") +
  geom_boxplot() +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  theme(strip.background = element_blank()) +
  geom_line(aes(group=paired, color = factor(smalldel)), linetype=2) +
  scale_color_manual(values = c("gray80", "red")) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_n_text()

pdf(file = "figures/spitzer_fig5b.pdf", height = 4, width = 4, useDingbats = FALSE)
small_del_plot
dev.off()

### ### ### ### ### ###
### Figure 5c - waterfall plot for hypoxia state with genetic annotations
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

malignant_freq_hq <- malignant_freq %>% 
  mutate(Patient =  sapply(strsplit(ID, "T"), "[[", 1),
         Timepoint= paste0("T", sapply(strsplit(ID, "T"), "[[", 2))) %>% 
  filter(Patient%in%longitudinal_malignant_cell_analysis$Patient, Timepoint!="T3")

small_deletion_values <- longitudinal_variant_del %>% 
  inner_join(linker_file_cases, by="case_barcode") %>% 
  filter(Patient%in%cosmic_sbs_v3_care_pairs$Patient) %>% 
  filter(!case_barcode%in%acquired_hypermutation$case_barcode)

sub_dat <- malignant_freq_hq %>% 
  filter(Patient%in%small_deletion_values$Patient,
         State%in%c("Hypoxia")) %>% 
  mutate(phenotype = ifelse(Patient%in%unique(c(sbs21_wide$Patient, small_del_cases_df$Patient)), "Small deletion/SBS21 increase", "Other"),
         smalldel = ifelse(Patient%in%unique(small_del_cases_df$Patient), "Small deletion", "Other"),
         sbs21 = ifelse(Patient%in%unique(sbs21_wide$Patient), "SBS21 increase", "Other")) %>% 
  dplyr::select(Patient, Timepoint, State, freq, phenotype, smalldel, sbs21) %>% 
  pivot_wider(names_from = Timepoint, values_from = freq) %>% 
  dplyr::select(Patient:sbs21, hypoxia_t1 = T1, hypoxia_t2 = T2) %>% 
  mutate(diff = (hypoxia_t2-hypoxia_t1)*100)


sub_dat <- sub_dat[order(sub_dat$diff, decreasing=TRUE),]
sub_dat <- sub_dat %>% 
  mutate(Patient = as_factor(Patient)) %>%
  mutate(Patient = fct_relevel(Patient, sub_dat$Patient))
sub_dat$fill <- sub_dat$diff > 0


gg_top <- ggplot(sub_dat, aes(x=Patient, y = diff, fill = factor(fill))) +
  geom_bar(stat="identity") +
  theme_classic() +
  scale_fill_manual(values=c("#27408B","#CD4F39")) +
  labs(title = "Hypoxia state", y = "Longitudinal change (%)") +
  theme(plot.title = element_text(size=14, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y= element_text(size=14),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")

small_deletion_values_filt <- small_deletion_values %>% 
  filter(Patient%in%sub_dat$Patient) %>% 
  mutate(Patient = as_factor(Patient)) 

small_deletion_values_filt$Patient <- factor(small_deletion_values_filt$Patient, levels= levels(sub_dat$Patient))


gg_middle <- ggplot(small_deletion_values_filt, aes(x=Patient, y=1, fill=mf_r)) +
  geom_tile(color = "white")+
  scale_fill_viridis(direction=-1, option="plasma") +
  labs(y = "",  fill="Small del\nT2", x="") + 
  plot_theme +
  theme(axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x= element_blank(),)

sbs21_wide_all <- cosmic_sbs_v3_care %>% 
  filter(Patient%in%cosmic_sbs_v3_care_pairs$Patient, SBS.Sig.max=="SBS21") %>% 
  dplyr::select(SBS.Sig.max, freq, Patient, Timepoint) %>% 
  pivot_wider(names_from = Timepoint, values_from = freq) %>% 
  mutate(sbs21_diff = T2-T1) %>% 
  filter(Patient%in%sub_dat$Patient)
sbs21_wide_all$Patient <- factor(sbs21_wide_all$Patient, levels= levels(sub_dat$Patient))

gg_bottom <- ggplot(sbs21_wide_all, aes(x=Patient, y=1, fill=T2)) +
  geom_tile(color = "white")+
  scale_fill_viridis(direction=-1, option="viridis") +
  labs(y = "SBS21", fill="SBS21\nT2") + 
  plot_theme +
  theme(axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 14, hjust = 1), text=element_text(size=14))

sub_dat$Patient <- factor(sub_dat$Patient, levels= levels(sub_dat$Patient))

gg_extra <- ggplot(sub_dat, aes(x=Patient, y = 1, fill = factor(smalldel))) +
  geom_tile() +
  labs(y="Small\ndel.") +
  plot_theme +
  scale_fill_manual(values=c("white","black"),na.value="#E5E5E5") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x= element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")

gg_extra2 <- ggplot(sub_dat, aes(x=Patient, y = 1, fill = factor(sbs21))) +
  geom_tile() +
  labs(y="SBS21") +
  plot_theme +
  scale_fill_manual(values=c("white","black"),na.value="#E5E5E5") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x= element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")

# Added a continuous variable for small deletion burden specific to T2.
egg::ggarrange(gg_top, gg_middle, gg_extra, gg_extra2, nrow = 4, heights = c(1.25, 0.25, 0.25, 0.25))

pdf("spitzer_fig5c.pdf", width = 6, height = 4.5, useDingbats = FALSE)
egg::ggarrange(gg_top, gg_extra, gg_extra2, nrow = 3, heights = c(1.5, 0.25, 0.25))
dev.off()

# Assign statistical significance to the longitudinal hypoxia changes in those samples with an acquired genetic alteration
wilcox_input <- malignant_freq_hq %>% 
  filter(Patient%in%small_deletion_values$Patient,
         State%in%c("Hypoxia")) %>% 
  mutate(phenotype = ifelse(Patient%in%unique(c(sbs21_wide$Patient, small_del_cases_df$Patient)), "Small deletion/SBS21 increase", "Other"),
         smalldel = ifelse(Patient%in%unique(small_del_cases_df$Patient), "Small deletion", "Other"),
         sbs21 = ifelse(Patient%in%unique(sbs21_wide$Patient), "SBS21 increase", "Other")) %>% 
  dplyr::select(Patient, Timepoint, State, freq, phenotype, smalldel, sbs21) 

ggplot(wilcox_input %>% 
         filter(smalldel=="Small deletion"), aes(x=Timepoint, y=freq*100)) +
  geom_boxplot() +
  geom_line(aes(group=Patient), color="gray70", linetype=2) +
  stat_compare_means(method = "wilcox", paired = TRUE) +
  plot_theme

ggplot(wilcox_input %>% 
         filter(sbs21=="SBS21 increase"), aes(x=Timepoint, y=freq*100)) +
  geom_boxplot() +
  geom_line(aes(group=Patient), color="gray70", linetype=2) +
  stat_compare_means(method = "wilcox", paired = TRUE) +
  plot_theme

# There's also a significant difference when comparing the longitudinal difference in Hypoxia between those samples with acquired small deletions versus all others.
ggplot(sub_dat, aes(x=smalldel, y=diff)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox") +
  plot_theme


### ### ### ### ### ###
### Figure 5d - GLASS ssGSEA for hypoxia signature in longitudinal pairs with acquired small deletions
### ### ### ### ### ###
# GLASS bulk RNA sequencing samples with both DNA and RNA that were classified as IDH-wildtype were scored for the Hypoxia (CARE) metaprogram
glass_dna_rna_smalldel <- read.delim("data/glass_longitudinal_hypoxia_ssgsea.txt", sep = "\t", header = TRUE)

pdf("figures/spitzer_fig5d_glass_hypoxia_scoring.pdf", width = 5, height = 4, useDingbats = FALSE)
ggplot(glass_dna_rna_smalldel, aes(x = timepoint, y =hypoxia_ss_gsea,  fill = timepoint)) + 
  geom_line(aes(group=paired), color="gray70", linetype=2) +
  geom_boxplot() +
  scale_fill_manual(values=c("T1"= "#66C2A5", "T2"="#FC8D62")) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.format", method = "wilcox", paired = TRUE) +
  labs(y="Bulk RNA ssGSEA\nsnRNA Hypoxia metaprogram", x = "GLASS - longitudinal pairs") +
  theme(strip.background = element_blank()) +
  stat_n_text() +
  facet_grid(.~type, scales="free")
dev.off()

### ### ### ### ### ###
### Figure 5e - transcriptional distance versus genetic distance
### ### ### ### ### ###
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
# Pearson correlation coefficient= 0.421 and p-value = 0.008465.
cor.test(genetic_transcript_dist_care$scaled_genetic_distance, genetic_transcript_dist_care$Dist01, method = "p")

pdf("figures/spitzer_fig5e_genetic_vs_transcriptomic_distance.pdf", width = 5, height = 4, useDingbats = FALSE)
ggplot(genetic_transcript_dist_care, aes(x=scaled_genetic_distance, y=Dist01)) +
  stat_cor(method = "pearson", cor.coef.name = "R") +
  geom_point() +
  plot_theme +
  labs(x = "Genetic distance T1-T2\n(SNV+CNA diff.)", y="Transcriptional distance T1-T2\n(Euclidean distance metric)") +
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "bottom")
dev.off()


### END ###