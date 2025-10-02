#################################
## Title: Extended Data Figure 2 plots (Spitzer et al. EDF2b, 2d, 2e, and 2f panels)
## Date: 2025.02.10
## Author: Kevin Johnson
## Description: Produce paired box plots for more granular cell states across time points and assorted Wang et al comparative analyses
#################################

# Necessary packages
library(tidyverse) # v1.3.1
library(ggpubr) # v0.4.0
library(EnvStats) # v2.7.0 
library(openxlsx) # v4.2.5.2
library(ggdist) # v3.2.1
# Custom minimalist plotting theme
source("R/plot_theme.R")

## Relevant input data.
# The cells analyzed post-filtering.
celltype_md <- readRDS("data/celltype_meta_data_2025_01_08.RDS")
# Following derivation of non-malignant metaprograms and removal of some lymphocytes due to ambiguous nature of some population and "other" (i.e., not classified cells).
# For each cell type, "_other" categories are where we could not confidently assign an expression metaprogram to a particular cell.
celltype_md_mp <- readRDS("data/mp_extended_cell_type_classification.RDS")

##########
### EDF 2b - selected more granular cell state abundance differences
##########
celltype_md_mp_filt <- celltype_md_mp %>% 
  dplyr::select(CellID, AssignedState)

# For this analysis, we wanted to restrict to "true" primaries. That is, where the T1 sample was from the patient's first surgery due to concerns that there might be differences when T1 was a recurrence.
true_primary_patients <- celltype_md_mp %>% 
  dplyr::select(Sample, Patient, Timepoint, PvsR, ID) %>% 
  distinct() %>% 
  filter(PvsR=="Primary")

# Perform a left_join so that we have both the full annotation and also more granular annotation.
all_cells_subclass = celltype_md %>%
  left_join(celltype_md_mp_filt, by="CellID") 

all_cells_subclass$AssignedState <- as.character(all_cells_subclass$AssignedState)
all_cells_subclass$CellTypeRefined <- as.character(all_cells_subclass$CellType)
all_cells_subclass$CellTypeRefined[!is.na(all_cells_subclass$AssignedState)] <- all_cells_subclass$AssignedState[!is.na(all_cells_subclass$AssignedState)]

# Calculating the cell type abundance with more granular annotations
nonmalignant_freq <- all_cells_subclass %>% 
  group_by(ID, CellTypeRefined) %>%
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(ID, CellTypeRefined,
           fill = list(counts = 0, freq = 0)) %>% 
  mutate(Population = "All") %>% 
  dplyr::select(ID, State = CellTypeRefined, counts, freq, Population) 


all_tme_freq_hq <- nonmalignant_freq %>% 
  separate(ID, into = c("Patient", "Timepoint"), sep = "T", remove = FALSE) %>% 
  mutate(Timepoint = paste0("T", Timepoint)) %>% 
  dplyr::select(Patient, State, Timepoint, freq) %>%
  pivot_wider(names_from = Timepoint, values_from = freq) %>% 
  filter(!is.na(T1), !is.na(T2)) %>% 
  dplyr::select(-T3) %>% 
  filter(Patient%in%true_primary_patients$Patient)

all_tme_freq_hq_long <- all_tme_freq_hq %>% 
  pivot_longer(cols=c(T1, T2), names_to = "Timepoint", values_to = "freq") %>% 
  arrange(Patient, State) %>% 
  mutate(paired = rep(1:(n()/2), each=2)) 

# 56 patients and 47 distinct states
wilcox_results <- all_tme_freq_hq_long %>%
  group_by(State) %>%
  summarise(
    p_value = wilcox.test(freq[Timepoint == "T1"], freq[Timepoint == "T2"], paired = TRUE)$p.value, # warning for exact p-values
    .groups = "drop"
  ) %>% 
  mutate(p_adjusted = p.adjust(p_value, method = "fdr"))

# There are five significant results after correction for multiple hypothesis testing.
# We are not reporting on the "Malignant" broad cell abundance difference at T2 (done in other figures) nor are we reporting on "Oligodendrocyte_other" because it's not clear what this group represents. Probably, a general oligodendrocyte increase.
wilcox_results_sig <- wilcox_results %>% 
  filter(p_adjusted < 0.05)

# We are not reporting on the Malignant cell difference nor are we reporting on "Oligodendrocyte_other".
selected_states <- c("Oligodendrocyte_mature", "Oligodendrocyte_mature_IFN",  "Astrocyte_reactive")

# Set factor levels based on relative cell abundance - high to low
all_tme_freq_hq_long_filt <- all_tme_freq_hq_long %>% 
  filter(State%in%selected_states)
all_tme_freq_hq_long_filt$State <- factor(all_tme_freq_hq_long_filt$State, levels=c( "Oligodendrocyte_mature", "Astrocyte_reactive", "Oligodendrocyte_mature_IFN"))


pdf("figures/spitzer_edf2b_granular_tme_abundance.pdf",width = 6, height = 5, useDingbats = FALSE)
ggplot(all_tme_freq_hq_long_filt  %>% mutate(State = recode(State, `Oligodendrocyte_mature`="Oligo.\nmature",
                                                            `Oligodendrocyte_mature_IFN` = "Oligo.\ninterferon",
                                                            `Astrocyte_reactive` = "Astro.\nreactive"),
                                             recode(Timepoint, `T1`="Primary",
                                                    `T2`="Recurrence")), aes(x = Timepoint, y = freq*100)) + 
  geom_line(aes(group=paired), color="gray70", linetype=2) +
  geom_boxplot(aes(fill=State)) +
  geom_point() +
  scale_linetype_manual(values="dashed") +
  facet_grid(.~State, scales="free") +
  plot_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("Oligo.\nmature"= "#B3DE69",
                             "Oligo.\ninterferon" = "#B3DE69",
                             "Astro.\nreactive"= "#BFBADA",
                             "Inh. Neuron\nMGE SST"= "#FFED6F")) +
  facet_grid(.~State,  scales="free") +
  stat_compare_means(method = "wilcox", 
                     paired = TRUE, 
                     size = 4, 
                     label="p.format") +
  labs(x="Time point", y="Relative state abundance (% of all cells)") +
  theme(strip.background = element_blank()) 
dev.off()

##########
# EDF 2d - Clinical data availability for CARE vs Wang et al
##########
# Wang et al (Diaz Laboratory) Nature Cancer 2022 - Downloaded Supplemental Table 1 and saved as a text file.
clin_md_snrna <- read.delim("data/diaz_clinical_snrna_metadata_20240106.txt", header = TRUE, sep = "\t")

# The total number of patients with clinical data in the Wang et al cohort - 34 patients
wang_patients <- length(unique(clin_md_snrna$Pair.))

# Filter to the number of Pair. (i.e., Patients) with non-missing overall survival data and surgical interval data.
pair_survival <- clin_md_snrna %>% 
  filter(!is.na(Overall.survival)) %>% 
  dplyr::select(Pair., Overall.survival) %>% 
  distinct()

pair_interval <- clin_md_snrna %>% 
  filter(!is.na(Elapsed.time.to.recurrence)) %>% 
  dplyr::select(Pair., Elapsed.time.to.recurrence) %>% 
  distinct()

# CARE clinical data - from Spitzer et al.
care_md <- readWorkbook("data/spitzer_supptable1.xlsx", colNames = TRUE)
# Fixing one instance where NA was incorrectly coded
care_md$MGMT.protomoter.methylation.status[care_md$MGMT.protomoter.methylation.status=="Na"] <- NA

# Restrict to the first time point since that's what we care about in terms of treatment because we evaluate biology before (T1) and after (T2) treatment. Treatment following T2 is less/not relevant for snRNA analyses.
care_md_patient <- care_md %>% 
  filter(Time.point=="T1")

# What's the breakdown for available MGMT information.
table(care_md_patient$MGMT.protomoter.methylation.status) # 42/59
sum(!is.na(care_md_patient$MGMT.protomoter.methylation.status))/dim(care_md_patient)[1]

# Create a data frame on the fly with the shared variables between the cohorts besides MGMT, which can be an important predictive variable for treatment response.
features <- c("Patient age", "Patient age",
             "Patient sex", "Patient sex",
             "Survival status", "Survival status",
             "Survival time", "Survival time",
             "Surgical time", "Surgical time",
             "MGMT status", "MGMT status")

clinical_df_compare <- data.frame(
  Cohort = c("CARE", "Wang", "CARE", "Wang", "CARE", "Wang", "CARE", "Wang", "CARE", "Wang", "CARE", "Wang"),
  ClinicalFeature = features,
  # Patient age and sex were reported across both cohorts
  PercentAvail = c(100, 100, 100, 100, 
                   # Survival status (alive or dead) was not listed for Wang et al
                   100, 0, 
                   # Survival time - no missingness for CARE
                   100, n_distinct(pair_survival$Pair.)/wang_patients*100, 
                   # Surgical interval - no missingness for CARE
                   100, n_distinct(pair_interval$Pair.)/wang_patients*100, 
                   # MGMT methylation status - not reported for Wang et al.
                   sum(!is.na(care_md_patient$MGMT.protomoter.methylation.status))/dim(care_md_patient)[1]*100, 0) 
)

pdf("figures/spitzer_edf2d_wang_care_avail_clinical.pdf", width = 5, height = 6, useDingbats = FALSE)
ggplot(clinical_df_compare, aes(x = Cohort, y = ClinicalFeature, fill=PercentAvail)) +
  geom_tile() +
  scale_fill_gradient2("Percent of data\navailable") +
  labs(x = "Dataset", y = "Clinical feature") +
  plot_theme +
  theme(legend.position = "top")
dev.off()

#########
# EDF 2e - nFeature RNA cohort comparisons
#########
# There are differences in the number of genes expressed across the different cell types. There are also differences in the cell type
# abundance between the two cohort. Reduce to only malignant cells.
diaz_md <- readRDS("data/diaz_cell_md_20240820.RDS")
diaz_md_malignant <- diaz_md %>% 
  filter(CellType=="Malignant")

care_cell_md <- readRDS("data/care_full_gene_count_metadata.RDS")
malignant_md <- readRDS("data/malignant_meta_data_2025_01_08.RDS")

care_md_malignant_compare <- care_cell_md %>% 
  filter(CellID%in%malignant_md$CellID) %>% 
  mutate(cohort = "CARE dataset") %>% 
  dplyr::select(CellID, cohort, nCount_RNA, nFeature_RNA) %>%
  # Applying the same filters used in Wang et al filtering
  filter(nFeature_RNA > 500)

diaz_md_malignant_compare <- diaz_md_malignant %>% 
  mutate(cohort = "Wang et al dataset") %>% 
  dplyr::select(CellID, cohort, nCount_RNA, nFeature_RNA)

cohort_df <- bind_rows(care_md_malignant_compare, diaz_md_malignant_compare) 

pdf("figures/spitzer_edf2e_wang_care_gene_complexity.pdf", width = 5, height = 5, useDingbats = FALSE)
ggplot(cohort_df, aes(x = cohort, y = nFeature_RNA)) + 
  ggdist::stat_halfeye(adjust = .5, width = .75, justification = -.2, .width = 0, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.color = NA) +
  coord_cartesian(xlim = c(1.2, NA)) +
  xlab("") +
  ylab("nFeature RNA") +
  scale_y_continuous(breaks = seq(0, 12000, 1000)) +
  stat_compare_means(method="wilcox") +
  ylim(0,10000) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

# Attempt to extract exact pvalue
wilcox_res <- wilcox.test(cohort_df$nFeature_RNA ~ cohort_df$cohort)
exact_p_value <- wilcox_res$p.value
print(format(exact_p_value, scientific = TRUE))

#########
# EDF 2f - Compositional group differences across time points in Wang et al
#########
# n = 50 samples, n = 25 pairs with compositional cluster annotation based on cell type abundance.
composition_cluster_summary <- readRDS("data/wang_compositional_clusters.RDS")

# p = 0.08 for all groups (e.g., LMF-TAM)
fisher.test(table(composition_cluster_summary$Stage, composition_cluster_summary$composition_cluster))

# p = 0.04 for HMF vs IMF vs LMF (restricting to broader groups)
fisher.test(table(composition_cluster_summary$Stage, composition_cluster_summary$comp_clust_groups))

# Enumerating the breakdown by Stage and compositional cluster
compclust_input <- composition_cluster_summary %>% 
  group_by(Stage, composition_cluster) %>% 
  summarise(counts = n()) %>%
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(Stage, composition_cluster,
           fill = list(counts = 0, freq = 0))

pdf("figures/spitzer_edf2f_wang_compositional_clusters.pdf", width = 6, height = 5, useDingbats = FALSE)
ggplot(compclust_input, aes(fill=Stage, x=composition_cluster, y=freq*100)) + 
  geom_bar(position="dodge", stat="identity") +
  plot_theme +
  scale_fill_manual(values=c("Primary" = "#66c2a5", "Recurrent" = "#fc8d62")) +
  labs(x="Composition group", y="Proportion (%)") +
  theme(legend.position = "top")
dev.off()  

### END ###