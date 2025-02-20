#################################
## Title: Extended Data Figure 8a-e in Nomura et al
## Date: 2025.02.15
## Author: Kevin Johnson
## Description: Produce heatmap for enrichment scores for tumor-associated myeloid cells (TAMs) and astrocytes 
#################################

# Necessary packages
library(tidyverse) # v1.3.1
library(viridis) # v0.6.5
library(scales) # v1.3.0
library(openxlsx) # v4.2.5.2
library(topGO) # v2.50.0
library(org.Hs.eg.db) # v3.16.0
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
# Custom minimalist plotting theme
source("R/plot_theme.R")

# Load in the malignant and non-malignant metaprograms from Nomura Supplementary Tables
malignant_mps <- readWorkbook("data/nomura_supptables.xlsx", sheet = 2, startRow = 5, colNames = TRUE)
nm_mps <- readWorkbook("data/nomura_supptables.xlsx", sheet = 4, startRow = 5, colNames = TRUE)
total_gene_set <- readRDS("data/care_total_gene_set.RDS")

### ### ### ###
### EDF 8a
### ### ### ###
# For each non-malignant metaprogram group, run a TopGO analysis - start with myeloid signatures
macrophage_sigs <- nm_mps %>% 
  as.data.frame() %>% 
  # Note: the metaprograms were defined as "Macrophage" but we later recategorized these as "TAM".
  dplyr::select(starts_with("Macrophage"))
mp_list <- as.list(macrophage_sigs)

# Function for enrichment of NMF metaprogram genes against covered background.
selFun = function(x) {
  ifelse(x==1, TRUE, FALSE)
}

goList <- list()
for(i in 1:length(mp_list)){
  cat("\r", i)
  mysig <- mp_list[[i]]
  mp_name <- names(mp_list)[i]
  all_genes <- total_gene_set
  
  gene_list_mp <- ifelse(all_genes%in%mysig, 1, 0)
  names(gene_list_mp) <- all_genes
  
  # Functional enrichment of metaprogram signature
  mpGOdata <- new("topGOdata",
                  ontology = "BP", # Biological process ontology
                  allGenes = gene_list_mp, # Gene list
                  geneSel = selFun, # Function to indicate which genes make up this list
                  annot=annFUN.org,
                  mapping = 'org.Hs.eg.db', # The annotation package for the human genome
                  ID = 'symbol', # We're using gene symbols
                  nodeSize = 10)
  
  # Fishers test
  resultFisher <- runTest(mpGOdata, algorithm = "classic", statistic = "fisher")
  fishRes <- GenTable(mpGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
  fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")
  
  fishRes[,"MP"] <- mp_name
  goList[[i]] <- fishRes
}

goRes <- do.call(rbind, goList)

# Inspect the top hits
top_30_q_values <- goRes %>%
  group_by(MP) %>%
  top_n(30, -q.value) %>% 
  filter(q.value<0.05)

top_1_q_values <- goRes %>%
  group_by(MP) %>%
  top_n(1, -q.value) %>% 
  filter(q.value<0.05)

# To plot one GO term per pathway, select a top hit that was used to assist cell state annotation. 
# The selected significant GO Terms were based on top hits and how the metaprogram aligned with publicly available gene signatures.
view(goRes)
terms_to_plot <- c("nervous system development", "endocytosis", "cellular response to interleukin-1", "cytoplasmic translation",
                   "response to hypoxia", "protein-lipid complex subunit organization", "inflammatory response", "response to type I interferon",
                   "cell cycle", 
                   #"nervous system development",
                   "protein phosphorylation", "MHC class II protein complex assembly",
                   #"endocytosis", 
                   "cell-cell adhesion", "response to unfolded protein", "MHC protein complex assembly")
goRes_plot <- goRes %>%
  filter(Term%in%terms_to_plot)

mp_order <- c("Macrophage_MP_1",
              "Macrophage_MP_2",
              "Macrophage_MP_3",
              "Macrophage_MP_4",
              "Macrophage_MP_5",
              "Macrophage_MP_6",
              "Macrophage_MP_7",
              "Macrophage_MP_8",
              "Macrophage_MP_9",
              "Macrophage_MP_10",
              "Macrophage_MP_11",
              "Macrophage_MP_12",
              "Macrophage_MP_13", 
              "Macrophage_MP_14",
              "Macrophage_MP_15",
              "Macrophage_MP_16")

goRes_plot$MP <- factor(goRes_plot$MP, levels=mp_order)
goRes_plot$Term <- as.factor(goRes_plot$Term)
goRes_plot$Term <- factor(goRes_plot$Term, levels=terms_to_plot)

pdf("figures/nomura_edf8a_myeloid_mp_go_enrichment.pdf", width = 6, height = 4.5)
ggplot(goRes_plot, aes(x = MP, y = Term, fill = -log10(q.value))) +
  geom_tile() +
  scale_fill_gradient2(low ="white", high = "#756bb1") +
  labs(y = "Enriched GO Term", x = "TAM MPs", fill="Significance\n-log10(adj.p)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1), text=element_text(size=8),
        legend.position = "top")
dev.off()



### ### ### ###
### EDF 8b
### ### ### ###
# Repeat this process for astrocyte metaprogram annotation
astrocyte_sigs <- nm_mps %>% 
  as.data.frame() %>% 
  dplyr::select(starts_with("Astrocyte"))
mp_list <- as.list(astrocyte_sigs)

goList <- list()
for(i in 1:length(mp_list)){
  cat("\r", i)
  mysig <- mp_list[[i]]
  mp_name <- names(mp_list)[i]
  all_genes <- total_gene_set
  
  gene_list_mp <- ifelse(all_genes%in%mysig, 1, 0)
  names(gene_list_mp) <- all_genes
  
  # Functional enrichment of metaprogram signature
  mpGOdata <- new("topGOdata",
                  ontology = "BP", # Biological process ontology
                  allGenes = gene_list_mp, # Gene list
                  geneSel = selFun, # Function to indicate which genes make up this list
                  annot=annFUN.org,
                  mapping = 'org.Hs.eg.db', # The annotation package for the human genome
                  ID = 'symbol', # We're using gene symbols
                  nodeSize = 10)
  
  # Fishers test
  resultFisher <- runTest(mpGOdata, algorithm = "classic", statistic = "fisher")
  fishRes <- GenTable(mpGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
  fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")
  
  fishRes[,"MP"] <- mp_name
  goList[[i]] <- fishRes
}

goRes <- do.call(rbind, goList)

# Inspect the top hits
top_30_q_values <- goRes %>%
  group_by(MP) %>%
  top_n(30, -q.value) %>% 
  filter(q.value<0.05)

top_1_q_values <- goRes %>%
  group_by(MP) %>%
  top_n(1, -q.value) %>% 
  filter(q.value<0.05)

# Manually curated GO terms to plot based on top hits. Needed to evaluate broader enrichment since many GO terms top hits are similar.
terms_to_plot <- c("neuron differentiation", "response to axon injury", "regulation of neurotransmitter levels", "phagocytosis",
                   "glutamate receptor signaling pathway", "response to wounding")

goRes_plot <- goRes %>%
  filter(Term%in%terms_to_plot)

goRes_plot$Term <- as.factor(goRes_plot$Term)
goRes_plot$Term <- factor(goRes_plot$Term, levels=terms_to_plot)

pdf("figures/nomura_edf8b_astrocyte_mp_go_enrichment.pdf", width = 6, height = 4.5)
ggplot(goRes_plot, aes(x = MP, y = Term, fill = -log10(q.value))) +
  geom_tile() +
  scale_fill_gradient2(low ="white", high = "#756bb1") +
  labs(y = "Enriched GO Term", x = "Astrocyte MPs", fill="Significance\n-log10(adj.p)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1), text=element_text(size=8),
        legend.position = "top")
dev.off()

### ### ### ###
### EDF 8c
### ### ### ###
# Individual myeloid cells scored for both CARE-defined gene expression metaprograms and publicly available myeloid scores using AddModuleScore.
# Gene signatures used:
# Muller et al 2017, PMID: 29262845. Supplementary Table 4.
# Hara et al 2021, PMID: 34087162. Supplementary Table 3.
# Venteicher et al 2017, PMID: 28360267. Supplementary Table 3.
# Atunes et al 2021, PMID: 33782623. Supplementary Table 2 and Supplementary Table 8 (recurrent)
# Miller et al 2023. PMID: 37961527 (Pre-print). Preliminary differentially expressed genes by author provided by author prior to pre-print - exploratory use only.
# Score for each gene signature for myeloid cells.
tam_mp_scores <- readRDS("data/myeloid_gene_set_scores.RDS")

# Convert to a matrix and assess the Pearson correlation coefficient
tam_mp_scores_mat <- tam_mp_scores %>% 
  dplyr::select(Macrophage_MP_1:Macrophage_MP_16, Muller2017_MG_Markers:Antunes2021_recurrent_hypoxic, Johnson2021_neutrophils:Johnson2021_dc) %>% 
  as.matrix()
tam_cormat <- round(cor(tam_mp_scores_mat), 2)

tam_cormat_mut_df <- data.frame(tam_cormat[1:16, 17:34])
tam_cormat_mut_df$WT_MP <- rownames(tam_cormat_mut_df)
tam_cormat_long <- tam_cormat_mut_df %>% 
  pivot_longer(cols= c(Muller2017_MG_Markers:Antunes2021_recurrent_hypoxic, Johnson2021_neutrophils:Johnson2021_dc),
               names_to = "MP",
               values_to = "Cor")


# Generate p-values
p_values_df <- data.frame(
  WT_MP = character(0),
  MP = character(0),
  p_value = numeric(0)
)

# Calculate p-values for each pair of variables
for (i in 1:(ncol(tam_cormat) - 1)) {
  for (j in (i + 1):ncol(tam_cormat)) {
    # Calculate the correlation and p-value
    cor_test_result <- cor.test(tam_mp_scores_mat[, i], tam_mp_scores_mat[, j])
    
    # Append the results to the p_values_df dataframe
    p_values_df <- rbind(p_values_df, data.frame(
      WT_MP = colnames(tam_cormat)[i],
      MP = colnames(tam_cormat)[j],
      p_value = cor_test_result$p.value
    ))
  }
}

tam_p_values_df_long <- p_values_df %>% 
  filter(WT_MP%in%tam_cormat_long$WT_MP, MP%in%tam_cormat_long$MP)

tam_cormat_long_comb <- tam_cormat_long %>% 
  inner_join(tam_p_values_df_long, by=c("WT_MP", "MP")) %>% 
  mutate(adj_p_value = p.adjust(p_value, "BH"))


mp_order <- c("Macrophage_MP_9",
              "Macrophage_MP_1",
              "Macrophage_MP_12",
              "Macrophage_MP_16",
              "Macrophage_MP_4",
              "Macrophage_MP_10",
              "Macrophage_MP_7",
              "Macrophage_MP_13",
              "Macrophage_MP_2",
              "Macrophage_MP_11",
              "Macrophage_MP_3",
              "Macrophage_MP_14",
              "Macrophage_MP_5", 
              "Macrophage_MP_15",
              "Macrophage_MP_6",
              "Macrophage_MP_8")

tm_mp_order <- c("Antunes2021_primary_prolif",
                 "Muller2017_MG_Markers",
                 "Hara2021_microglia",
                 "Venteicher2017_microglia",
                 "Antunes2021_primary_tam2",
                 "Johnson2021_neutrophils",
                 "Muller2017_Mac_Markers",
                 "Hara2021_macrophage",
                 "Hara2021_MESlike",
                 "Venteicher2017_macrophage",
                 "Antunes2021_primary_tam1",
                 "Antunes2021_recurrent_ifn",
                 "Antunes2021_recurrent_transitory", 
                 "Antunes2021_recurrent_lipid", 
                 "Antunes2021_recurrent_hypoxic",
                 "Antunes2021_primary_monocytes",
                 "Antunes2021_primary_dc",
                 "Johnson2021_dc")


tam_cormat_long_comb$WT_MP <- factor(tam_cormat_long_comb$WT_MP, levels=mp_order)
tam_cormat_long_comb$MP <- factor(tam_cormat_long_comb$MP, levels=rev(tm_mp_order))

pdf("figures/nomura_edf8c_myeloid_gene_set_scores.pdf", width = 6, height = 4.5)
ggplot(data = tam_cormat_long_comb, aes(y=MP,x=WT_MP, fill = Cor)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1), text=element_text(size=8),
        legend.position = "top") +
  labs(y= "Public myeloid signatures", x= "TAM MPs") 
dev.off()


### ### ### ###
### EDF 8d
### ### ### ###
# Individual astrocytes scored for both CARE-defined gene expression metaprograms and publicly available myeloid scores using AddModuleScore.
# Sadick et al 2022, PMID: 35381189. 
# Liu et al 2023, PMID: 36931245. 
# Score for each gene signature for astrocytes
astrocyte_mp_scores <- readRDS("data/astrocyte_gene_set_scores.RDS")

astro_mp_scores_mat <- astrocyte_mp_scores %>% 
  dplyr::select(Astrocyte_MP_1:Astrocyte_MP_7, Sadick2022_C0_mature_astro = `C0_protective-sadick2022`, Sadick2022_C3_reactive_astro  = `C3_reactive-sadick2022`, Liu2023_fetal_astro = `Astrocyte-liu2023`, Liu2023_outer_radial_glia = `Outer.radial.glia-liu2023`)%>% 
  as.matrix()

astro_cormat <- round(cor(astro_mp_scores_mat), 2)

cormat_mut_df <- data.frame(astro_cormat[1:7, 8:11])
cormat_mut_df$WT_MP <- rownames(cormat_mut_df)
cormat_long <- cormat_mut_df %>% 
  pivot_longer(cols= c(Sadick2022_C0_mature_astro:Liu2023_outer_radial_glia),
               names_to = "MP",
               values_to = "Cor")


# Generate p-values
astro_p_values_df <- data.frame(
  WT_MP = character(0),
  MP = character(0),
  p_value = numeric(0)
)

# Calculate p-values for each pair of variables
for (i in 1:(ncol(astro_cormat) - 1)) {
  for (j in (i + 1):ncol(astro_cormat)) {
    # Calculate the correlation and p-value
    cor_test_result <- cor.test(astro_mp_scores_mat[, i], astro_mp_scores_mat[, j])
    
    # Append the results to the astro_p_values_df dataframe
    astro_p_values_df <- rbind(astro_p_values_df, data.frame(
      WT_MP = colnames(astro_cormat)[i],
      MP = colnames(astro_cormat)[j],
      p_value = cor_test_result$p.value
    ))
  }
}

astro_p_values_df$MP <- gsub("-", "\\.", astro_p_values_df$MP)
astro_p_values_df_long <- astro_p_values_df %>% 
  filter(WT_MP%in%cormat_long$WT_MP, MP%in%cormat_long$MP)

cormat_long_comb <- cormat_long %>% 
  inner_join(astro_p_values_df_long, by=c("WT_MP", "MP")) %>% 
  mutate(adj_p_value = p.adjust(p_value, "BH"))

pdf("figures/nomura_edf8d_astrocyte_gene_set_scores.pdf", width = 6, height = 4.5)
ggplot(data = cormat_long_comb, aes(y=MP,x=WT_MP, fill = Cor)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1), text=element_text(size=8),
        legend.position = "top") +
  labs(y= "Public astrocyte signatures", x= "Astrocyte MPs") 
dev.off()


### ### ### ###
### EDF 8e
### ### ### ###
# Calculate overlaps for each pair of columns using apply.
intersect_matrix <- apply(nm_mps, 2, function(x) {
  sapply(malignant_mps, function(y) length(intersect(x, y)))
})

# Reshape the data and specify order to be presented. We previously called TAMs "Macrophages", but revised it to be tumor-associated myeloid (TAMs).
intersect_matrix_melt <- reshape2::melt(intersect_matrix) %>% 
  mutate(CellType = sapply(strsplit(as.character(Var2), "_MP"), "[[", 1)) %>% 
  mutate(CellType = recode(CellType, `Macrophage` = "TAM"))

cell_type_order <- c("Astrocyte", "Oligodendrocyte", "OPC", "ExN", "InN", "Endothel", "Pericyte", "TAM", "Tcell", "Bcell")
intersect_matrix_melt$CellType <- factor(intersect_matrix_melt$CellType, levels=cell_type_order)

# Modify the malignant metaprograms to have a cleaner aesthetic.
intersect_matrix_melt$Var1 <- gsub("MP_", "MP", intersect_matrix_melt$Var1)
intersect_matrix_melt$Var1 <- gsub("_", " ", intersect_matrix_melt$Var1)
intersect_matrix_melt$Var1 <- gsub("^(\\S+)\\s+(\\S+)$", "\\1 (\\2)", intersect_matrix_melt$Var1)

# Specify the order of the x-axis.
intersect_matrix_melt$Var1 <- factor(intersect_matrix_melt$Var1, levels=c("MP1 (RP)", "MP2 (OPC)", "MP3 (CC)", "MP4 (AC)",
                                                                          "MP5 (Hypoxia)", "MP6 (MES)", "MP7 (NPC)", "MP8 (GPC)", "MP9 (ExN)",
                                                                          "MP10 (Stress1)", "MP11 (MIC)", "MP12 (LQ)", "MP13 (Cilia)",
                                                                          "MP14 (NRGN)", "MP15 (Stress2)"))


pdf(file = "figures/nomura_edf8e_tme_mp_malignant_mp_similarity.pdf", height = 7, width = 8, useDingbats = FALSE)
ggplot(data = intersect_matrix_melt %>% 
         filter(CellType%in%c("Astrocyte", "TAM")), aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
  facet_grid(CellType~.,  scales = "free_y", space="free") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top") +
  labs(x="Malignant MPs", y="Non-malignant MPs")
dev.off()

# Cross-check annotation.
celltype_md_mp <- readRDS("data/mp_extended_cell_type_classification.RDS")
celltype_md_mp_filt <- celltype_md_mp %>% 
  filter(CellType%in%c("TAM", "Astrocyte")) %>% 
  dplyr::select(CellType, AssignedMP, AssignedState, CellTypeRefined) %>% 
  distinct()
celltype_md_mp_filt


### END ###
