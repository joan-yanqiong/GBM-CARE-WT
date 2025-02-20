#################################
## Title: Wang et al T1 vs T2 malignant cell abundance and malignant cluster (Spitzer et al. EDF3b-c panels)
## Date: 2025.02.10
## Author: Kevin Johnson
## Description: Produce paired box plots for Wang et al malignant cell state abundance changes across time points and bar plots for malignant clusters
#################################

# Necessary packages
library(tidyverse) # v1.3.1
library(ggpubr) # v0.4.0
library(EnvStats) # v2.7.0 
# Custom minimalist plotting theme
source("R/plot_theme.R")

# Per tumor sample, malignant cell state abundance count and frequency plus additional metadata
diaz_malignant_md_out <- readRDS("data/EDF3b_wang_et_al_small.RDS")

# Set the order for how you want the cell states to appear
cell_state_order <- c("Unresolved", 
                      "Neuron-like",
                      "NPC-like",
                      "OPC-like",
                      "GPC-like",
                      "Hypoxia",
                      "MES-like",
                      "Stress",
                      "AC-like" ,
                      "Cilia-like")
diaz_malignant_md_out$State <- factor(diaz_malignant_md_out$State, levels=rev(cell_state_order))

pdf("figures/spitzer_edf3b_wang_et_al.pdf", width = 4.5, height = 4.5, useDingbats = FALSE)
ggplot(diaz_malignant_md_out, aes(x = Stage, y = freq*100)) + 
  geom_line(aes(group=UCSF_Pair), color="gray70", linetype=2) +
  geom_boxplot(aes(fill=State)) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 3, label="p.format") +
  labs(x="Wang et al", y="Relative malignant cell abundance (%)") +
  scale_fill_manual(values=c("Cilia-like" = "#AA2756", 
                             "AC-like" ="#DA5361", 
                             "MES-like"="#F77D58",
                             "Hypoxia"="gray30",
                             "Stress"="#FCB672",
                             "GPC-like"="#FFFFC7",
                             "OPC-like"="#E8F5A3",
                             "NPC-like"="#B2E1AB",
                             "Neuron-like"= "#78C8AF",
                             "Unresolved" = "gray70")) +
  facet_grid(.~State, scales="free") +
  theme(strip.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_n_text()
dev.off()

### ### ### ### ### ### ###
## EDF3c in Spitzer et al
### ### ### ### ### ### ###
# Read in the malignant group assignment data for Wang et al
malignant_cluster_summary <- readRDS("data/wang_malclust_input.RDS")

malclust_input <- malignant_cluster_summary %>% 
  group_by(Stage, malignant_cluster) %>% 
  summarise(counts = n()) %>%
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(Stage, malignant_cluster,
           fill = list(counts = 0, freq = 0))

pdf("figures/spitzer_edf3c_wang_malignant_clusters.pdf", width = 6, height = 5, useDingbats = FALSE)
ggplot(malclust_input, aes(fill=Stage, x=malignant_cluster, y=freq*100)) + 
  geom_bar(position="dodge", stat="identity") +
  plot_theme +
  scale_fill_manual(values=c("Primary" = "#66c2a5", "Recurrent" = "#fc8d62")) +
  labs(x="Malignant group", y="Proportion (%)", fill="Wang et al") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()  


### END ###