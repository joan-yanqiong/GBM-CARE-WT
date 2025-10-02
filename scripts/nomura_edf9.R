#################################
## Title: Extended Data Figure 9 in Nomura et al
## Date: 2025.02.15
## Author: Kevin Johnson
## Description: Produce a large correlation heatmap for the cellular state abundances across malignant and non-malignant compartments
#################################

# Necessary packages
library(tidyverse) # v1.3.1
library(magrittr)  # v2.0.3
library(corrplot) # v0.92
library(Hmisc) # v4.7.0

# Load the full cell annotations across broad and more granular classifications
malignant_md <- readRDS("data/malignant_meta_data_2025_01_08.RDS")
celltype_md <- readRDS("data/celltype_meta_data_2025_01_08.RDS")
celltype_md_mp <- readRDS("data/mp_extended_cell_type_classification.RDS")

# Note that there are no malignant cells for one sample - P58T3
celltype_md_P58T3 <- celltype_md %>% 
  filter(Sample=="NL022")
table(celltype_md_P58T3$CellType)

### ### ### ### ### ###
### Combine malignant and non-malignant frequencies
### ### ### ### ### ###
# Examine malignant frequency of cycling cells - n = 120.
malignant_freq_cc <- malignant_md %>% 
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

malignant_freq = malignant_md %>%
  group_by(ID, State) %>%
  mutate(State = paste0(State, "-like")) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(ID, State,
           fill = list(counts = 0, freq = 0)) %>%
  mutate(Population = "malignant") %>% 
  # I am adding Cycling malignant cells as a feature in this matrix.
  bind_rows(malignant_freq_cc) %>% 
  # Recoding states to be consistent with the rest of the manuscript
  mutate(State = recode(State, `Hypoxia-like` = "Hypoxia",
                        `Stress-like` = "Stress",
                        `Neuron-like` = "NEU-like"))

# Note that malignant and non-malignant cell frequencies are calculated separately to adjust for potential differences in tumor purity.
nonmalignant_freq = celltype_md_mp %>%
  group_by(ID, AssignedState) %>%
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(ID, AssignedState,
           fill = list(counts = 0, freq = 0)) %>% 
  mutate(Population = "TME") %>% 
  dplyr::select(ID, State = AssignedState, counts, freq, Population)

# This object represents three different proportions collapsed into 1 df.
all_cell_freq <- malignant_freq %>% 
  bind_rows(nonmalignant_freq) 

# Remove any sample with less than 50 tumor cells
malignant_cell_number <- malignant_md %>% 
  group_by(ID) %>% 
  summarise(counts = n()) 

# Set threshold for samples that should be kept in an analysis at 50 malignant cells. This threshold is used as a rough approximation to reduce noise in the analysis.
malignant_cell_number_keep <- malignant_cell_number %>% 
  filter(counts >50)

# Filter to a high-quality set.
all_cell_freq_hq <- all_cell_freq %>% 
  filter(ID%in%malignant_cell_number_keep$ID)

# Create a structure that can be fed to corrplot.
state_freq <- all_cell_freq_hq %>% 
  dplyr::select(ID, State, freq) %>% 
  pivot_wider(names_from = State, values_from = freq)

# Compute the correlation matrix
dim(state_freq)
state_freq_cor <- cor(as.matrix(state_freq[ ,2:56]))
state_freq_cor_pval <- rcorr(as.matrix(state_freq[ ,2:56]))

# Plot with all correlations shown.
corrplot(state_freq_cor, method = 'circle', type = 'lower',
         tl.col = "black",
         tl.cex = 0.5,
         order = 'FPC', # 'FPC' for the first principal component order.
         insig='blank', number.cex = 0.8, diag=FALSE, col=rev(COL2("RdBu")))

# Include a p-value matrix in order to restrict to showing only significant associations. All else will be blank.
pdf(file = "figures/nomura_edf9_state_correlation_heatmap_large.pdf", height = 6, width = 6, useDingbats = FALSE)
corrplot(state_freq_cor, p.mat = state_freq_cor_pval$P, method = 'circle', type = 'lower',
         tl.col = "black",
         tl.cex = 0.5,
         order = 'FPC', # 'FPC' for the first principal component order.
         insig='blank', number.cex = 0.8, diag=FALSE, col=rev(COL2("RdBu")))
dev.off()

### END ###