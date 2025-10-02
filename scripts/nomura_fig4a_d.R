#################################
## Title: Figure 4 panels a-d in Nomura et al
## Date: 2025.02.15
## Author: Kevin Johnson
## Description: Examination of more granular tumor microenvironmental cell state abundances with baseline profiles and malignant cell state abundance.
#################################

# Necessary packages
library(tidyverse) # v1.3.1
library(magrittr)  # v2.0.3
library(corrplot) # v0.92
library(Hmisc) # v4.7.0
library(igraph) # v1.3.1
library(openxlsx) # v4.2.5.2
library(ggpubr) # v0.4.0
library(scales) # v1.3.0 
source("R/plot_theme.R")

# Load the full cell annotations across broad and more granular classifications
malignant_md <- readRDS("data/malignant_meta_data_2025_01_08.RDS")
celltype_md <- readRDS("data/celltype_meta_data_2025_01_08.RDS")
celltype_md_mp <- readRDS("data/mp_extended_cell_type_classification.RDS")
nm_mps <- readWorkbook("data/nomura_supptables.xlsx", sheet = 4, startRow = 5, colNames = TRUE)

### ### ### ### ### ###
# Figure 4a - donut chart for cell number
## ### ### ### ### ###
# Extract across cohort cell number values that were used as input for deriving gene expression metaprograms
table(celltype_md$CellType)
table(celltype_md_mp$AssignedState)
oligo_num <- as.numeric(table(celltype_md$CellType)[[2]])
astro_num <- as.numeric(table(celltype_md$CellType)[[4]])
opc_num <- as.numeric(table(celltype_md$CellType)[[10]])
exn_num <- as.numeric(table(celltype_md$CellType)[[5]])
inh_num <- as.numeric(table(celltype_md$CellType)[[6]])
endo_num <- as.numeric(table(celltype_md$CellType)[[7]])
pericyte_num <- as.numeric(table(celltype_md$CellType)[[8]])
tam_num <- as.numeric(table(celltype_md$CellType)[[3]])
tcell_num <- as.numeric(table(celltype_md$CellType)[[9]])-as.numeric(table(celltype_md_mp$AssignedState)[[37]])
bcell_num <- as.numeric(table(celltype_md_mp$AssignedState)[[37]])

# Create non-malignant cell proportion breakdown:
data <- data.frame(
  category=c("Oligodendrocyte", "Astrocyte", "OPC", "Excitatory neuron", "Inhibitory neuron", "Endothelial", "Pericyte", "TAM", "Lymphocyte\nTcell", "Lymphocyte\nBcell"),
  count=c(oligo_num, astro_num, opc_num, exn_num, inh_num, endo_num, pericyte_num, tam_num, tcell_num, bcell_num)
)

colnames(nm_mps)
# Manually enter the number of gene expression metaprograms detected in CARE.
num_mps <- c(9, 7, 1, 6, 4, 6, 5, 16, 6, 2)
names(num_mps) <- num_mps
data_mps <- cbind(data, num_mps)
# The number of metaprograms detected is positively correlated with cell number as expected. Greater power to detect expression variability with more cells.
cor.test(data_mps$count, data_mps$num_mps)
ggplot(data_mps, aes(x=count, y=num_mps)) + 
  geom_point()

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages 
data$ymax <- cumsum(data$fraction)


data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Provide an informative label
data$label <- paste0(data$category)

pdf(file = "figures/nomura_fig4a_tme_mp_donut_chart.pdf", height = 5, width = 5, useDingbats = FALSE)
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_manual(values=c("Oligodendrocyte" = "#B3DE69", "TAM" = "#80B1D3", "Astrocyte" = "purple", "Excitatory neuron" = "#BC80BD", 
                             "Inhibitory neuron" = "#FFED6F", "Endothel" = "#FCCDE5", "Pericyte" = "#FFFFB3", 
                             "Lymphocyte\nTcell" = "#8DD3C7", "Lymphocyte\nBcell" = "#BEBADA", "OPC" = "#FDB462")) +
  coord_polar(theta = "y", start = 0)+
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")
dev.off()

### ### ### ### ### ###
### Figure 4b - Selected TME abundance values across different baseline profiles (BPs)
### ### ### ### ### ###
# Load the state-controlled baseline profiles
scp <- readWorkbook("data/nomura_supptables.xlsx", sheet = 5, startRow = 4, colNames = TRUE)
# Remove the one sample without values due to missing malignant cells
scp_slim <- scp %>% 
  dplyr::select(SampleID, BP) %>% 
  filter(!is.na(BP))

celltype_md_mp_tme_clean = celltype_md_mp %>% 
  mutate(CellType = recode(CellType, `Endothel` = "Endothelial"))

nm_celltype_md_mp <- celltype_md_mp %>% 
  mutate(CellType = recode(CellType, `Endothel` = "Endothelial")) %>% 
  filter(CellType!="Malignant")

nm_scp <- nm_celltype_md_mp %>%
  group_by(ID, AssignedState) %>%
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(ID, AssignedState,
           fill = list(counts = 0, freq = 0)) %>% 
  inner_join(scp_slim, by=c("ID"="SampleID")) %>%
  # Removing "Mixed" since this is poorly defined.
  filter(BP!="Mixed")
# 97 samples left in this analysis
n_distinct(nm_scp$ID)
sum(unique(nm_scp$ID)%in%unique(malignant_md_counts$ID))

pdf("figures/nomura_fig4b_select_tme_states_bp.pdf", width = 7, height = 4.5, useDingbats = FALSE)
ggplot(nm_scp %>% 
         filter(AssignedState%in%c("Astrocyte_mature", "TAM_MD", "Endothelial_inflammatory", "TIL_CD4_Tcells")), aes(x=BP,  y=freq*100)) +
  geom_boxplot() +
  geom_point() +
  labs(y = "Sample relative\nTME abundance (%)", fill="Cell state", x = "") +
  stat_compare_means(method="kruskal", label = "p.format") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=10)) +
  facet_grid(.~AssignedState, scales = "free", space = "free") 
dev.off()


### ### ### ### ### ###
### Figure 4c - network graph to visualize the cell abundance correlations 
### ### ### ### ### ###
# Combine malignant and non-malignant frequencies
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
  bind_rows(malignant_freq_cc) %>% 
  mutate(State = recode(State, `Hypoxia-like` = "Hypoxia",
                        `Stress-like` = "Stress",
                        `Neuron-like` = "NEU-like"))


nonmalignant_freq = celltype_md_mp %>%
  group_by(ID, CellTypeRefined) %>%
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(ID, CellTypeRefined,
           fill = list(counts = 0, freq = 0)) %>% 
  mutate(Population = "TME") %>% 
  dplyr::select(ID, State = CellTypeRefined, counts, freq, Population) 

## This object represents three different proportions collapsed into 1 df.
all_cell_freq <- malignant_freq %>% 
  bind_rows(nonmalignant_freq) 

## Remove any sample with less than 50 tumor cells
malignant_cell_number <- malignant_md %>% 
  group_by(ID) %>% 
  summarise(counts = n()) 

## Set threshold for samples that should be kept in an analysis at 50 malignant cells.
malignant_cell_number_keep <- malignant_cell_number %>% 
  filter(counts >50)

all_cell_freq_hq <- all_cell_freq %>% 
  filter(ID%in%malignant_cell_number_keep$ID) 

# Create a corrplot:
state_freq <- all_cell_freq_hq %>% 
  dplyr::select(ID, State, freq) %>% 
  pivot_wider(names_from = State, values_from = freq)

# Assess the p-value and correlation cutoff
dim(state_freq)
cor_matrix <- cor(as.matrix(state_freq[ ,2:25]))

# Calculate p-values for the correlations using cor.mtest
p_values <- cor.mtest(as.matrix(state_freq[ ,2:25]))$p

# Set the significance threshold (e.g., 0.05) so that the graph is not overloaded
significance_threshold <- 0.05

# Initialize variables to store the minimum positive correlation and its indices
min_positive_corr <- Inf
row_index <- 0
col_index <- 0

# Iterate through the correlation matrix and p-values
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    if (i != j && p_values[i, j] <= significance_threshold && cor_matrix[i, j] > 0) {
      # If the correlation is positive and significant, update the minimum if needed
      if (cor_matrix[i, j] < min_positive_corr) {
        min_positive_corr <- cor_matrix[i, j]
        row_index <- i
        col_index <- j
      }
    }
  }
}

# Check if a significant positive correlation was found
if (min_positive_corr == Inf) {
  cat("No significant positive correlation below 0.001 found.\n")
} else {
  cat("Minimum positive correlation below 0.001:", min_positive_corr, "\n")
  cat("Corresponding row index:", row_index, "\n")
  cat("Corresponding column index:", col_index, "\n")
}


# Keep only high correlations
state_freq_cor <- cor(as.matrix(state_freq[ ,2:25]))
state_freq_cor_pval <- rcorr(as.matrix(state_freq[ ,2:25]))

# find cells with p-values > 0.05 and replace corresponding
# correlations coefficients with zero
state_freq_cor[state_freq_cor_pval$P > 0.0005] <- 0
state_freq_cor[state_freq_cor<0] <- 0

network <- graph_from_adjacency_matrix(state_freq_cor, weighted=T, mode="undirected", diag=F)
plot(network)

# Remove any nodes that are not connected  
Isolated = which(degree(network)==0)
network2 = delete.vertices(network, Isolated)

# Set the node colors in the graph based on the defined colors
V(network2) 
node_colors2 <- c(rep("#FB8072", 10), "#FCCDE5" , "#BC80BD", "#FFED6F", "#8DD3C7", "#FFFFB3", "#CDE6C4", rep("#80B1D3", 5))
V(network2)$color <- node_colors2

# Edge colors to represent correlation strength
color_palette <- colorRampPalette(c("gray80", "gray50", "gray20"))
edge_colors <- color_palette(10)[cut(E(network2)$weight, breaks = 10)]

edge_width <- c(3,6,9)[cut(E(network2)$weight, breaks = 3)]

# Plot the graph with node colors
set.seed(43)
pdf(file = "figures/nomura_fig4c_n111_samples_more_than_50_malignant_cells.pdf", height = 8, width = 8, useDingbats = FALSE)
plot(network2,
     vertex.label.cex=1,
     vertex.label.font=1,
     vertex.frame.color=NA,
     vertex.label.color="black",
     edge.color = "gray50",
     edge.width=edge_width)
dev.off()

### ### ### ### ### ###
### Figure 4d - Correlation heatmap showing more granular cell states and associations with malignant state abundance
### ### ### ### ### ###
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

# Produce a focused heatmap on hypoxia - extracting key features to visualize the association across different cell type categories.
colnames(state_freq)[c(c(2:10), 15, 16, 18, 42, 43)]
cormat <- round(cor(state_freq[ ,c(c(2:10), 15, 16, 18, 42, 43)]), 2)
cormat_mut_df <- data.frame(cormat[1:9, 10:14])

cormat_mut_df$State <- rownames(cormat_mut_df)
cormat_long <- cormat_mut_df %>% 
  pivot_longer(cols= c(Astrocyte_reactive:TAM_inflammatory_UPR),
               names_to = "MP",
               values_to = "Cor")

state_order <- c("AC-like",
                 "Cilia-like",
                 "GPC-like",
                 "OPC-like",
                 "NPC-like",
                 "NEU-like",
                 "MES-like",
                 "Stress",
                 "Hypoxia")
cormat_long$State <- factor(cormat_long$State, levels=rev(state_order))

pdf("figures/nomura_fig4d_tme_hypoxia_heatmap.pdf", width = 6, height = 4.5, useDingbats = FALSE)
ggplot(data = cormat_long, aes(x=MP,y=State, fill = Cor)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.65, 0.65), oob=squish, space = "Lab", 
                       name="Pearson Correlation") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=10)) +
  labs(y= "GBM cell states", x= "TME cell states") +
  theme(legend.position="top")
dev.off()


### END ###