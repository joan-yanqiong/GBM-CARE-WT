#################################
## Title: Figure 3 in Nomura et al - run the baseline profile analysis
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: This script generates the state-controlled baseline profiles that characterize the inter-tumor heterogeneity based on genes that
## recur across state-specific pseudo-bulk profiles
#################################

####################################################################################################################################
# Generate the PB profiles
####################################################################################################################################

d <- mdata
length(table(d$Sample))

# Generate profiles from the 7 principal cellular states
d_state_stats <- d %>%
  filter(State %in% c("AC", "MES", "Hypoxia", "GPC", "OPC", "NPC", "Neuron")) %>%
  group_by(Patient, Timepoint, Sample, State) %>%
  summarise(n = n())

n_cells <- 25

states <- c("AC", "MES", "Hypoxia", "GPC", "OPC", "NPC", "Neuron")

d_state_stats <- d_state_stats %>%
  filter(n >= n_cells)

# Remove "junk genes" (i.e. pseudo-genes, antisense genes etc)
junk_genes <- c(rownames(umi_data_all[[1]])[grep("\\.", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("-AS*", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("LINC", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("^RP[S|L]", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("^MT-", rownames(umi_data_all[[1]]))])

valid_genes <- rownames(umi_data_all[[1]])[rownames(umi_data_all[[1]]) %ni% junk_genes]

# Generate the pseudo-bulk profiles
state_pb_profiles <- lapply(unique(d_state_stats$State), function(st) {
  
  print(st)
  
  d_state_stats_tmp <- d_state_stats %>%
    filter(State == st)
  
  res <- lapply(d_state_stats_tmp$Sample, function(sname) {
    
    m <- umi_data_all[[sname]]
    
    d_tmp <- d %>%
      filter(Sample == sname, State == st)
    
    m <- m[valid_genes, d_tmp$CellID]
    
    m <- umi2upm(m)
    
    m <- log2(rowMeans(m) + 1)
    
    return(m)
  })
  res <- do.call(cbind, res)
  colnames(res) <- paste0(d_state_stats_tmp$Sample, "_", st)
  
  return(res)
})
names(state_pb_profiles) <- unique(d_state_stats$State)

####################################################################################################################################
# Select genes for the analysis
####################################################################################################################################

complexity <- apply(do.call(cbind, state_pb_profiles), 2, function(x) sum(x > 0))
n_cells <- setNames(d_state_stats$n, paste0(d_state_stats$Sample, "_", d_state_stats$State))

# Remove genes that don't have an entrezid
valid_genes <- tibble(Gene = rownames(state_pb_profiles[[1]]),
                      Entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, Gene, 'ENTREZID', 'SYMBOL'))
valid_genes <- valid_genes %>%
  filter(!is.na(Entrez))

rm_list <- lapply(state_pb_profiles, rowMeans)
rm_list <- lapply(rm_list, function(x) x[names(x) %in% valid_genes$Gene])

rm_list <- lapply(rm_list, function(x) x[x > 4])

rm <- rowMeans(do.call(cbind, state_pb_profiles))
rm <- rm[names(rm) %in% valid_genes$Gene]

# Variance is computed within each state to capture the genes that change within state (rather than across states)
rv <- lapply(state_pb_profiles, function(x) apply(x, 1, var))
rv <- do.call(cbind, rv)
rv <- apply(rv, 1, median)
rv <- rv[names(rv) %in% valid_genes$Gene]
rv <- rv[names(rm)]

sg <- tibble(Gene = names(rv), V = rv, M = rm)

# Select genes with mean expression > 1 abd variance > 2.5 (across all pseudo-bulk profiles)
sg <- sg %>%
  filter(M > 1, V > 2.5)
dim(sg)

sg <- sg$Gene

####################################################################################################################################
# Compute PCA per state to detect the genes contributing most variance within each state
####################################################################################################################################

pca_res_list <- lapply(names(state_pb_profiles), function(st) {
  
  print(st)
  
  pbm <- state_pb_profiles[[st]][sg, ]
  dim(pbm)
  
  # Center within state  
  pbm <- rowcenter(pbm)
  
  pca_res <- pcaMethods::pca(object = t(pbm), nPcs = 10, scale = "none", center = F)
  
  return(pca_res)
})
names(pca_res_list) <- names(state_pb_profiles)

pca_scores_tbl <- do.call(rbind,
                          lapply(pca_res_list, function(x) as_tibble(pcaMethods::scores(x), rownames = "ID")))
pca_scores_tbl$Sample <- sapply(strsplit(pca_scores_tbl$ID, "_"), function(x) x[[1]][1])
pca_scores_tbl$State <- sapply(strsplit(pca_scores_tbl$ID, "_"), function(x) x[[2]][1])

####################################################################################################################################
# Generate two signatures from each (State, PC)
####################################################################################################################################

npcs <- 3
n_genes <- 50

pca_loadings_tbl <- lapply(names(pca_res_list), function(st) {
  
  print(st)
  
  pca_res <- pca_res_list[[st]]
  
  lds <- pcaMethods::loadings(pca_res)
  dim(lds)
  
  lds_res <- lapply(1:npcs, function(i) {
    res <- tibble(V1 = lds[, i] %>% sort(decreasing = T) %>% head(n_genes) %>% names(),
                  V2 = lds[, i] %>% sort(decreasing = F) %>% head(n_genes) %>% names())
    
    colnames(res) <- paste0(st, "_PC", i, c("_HIGH", "_LOW"))
    
    return(res)
  })
  lds_res <- do.call(cbind, lds_res)
  
  return(lds_res)
})
pca_loadings_tbl <- do.call(cbind, pca_loadings_tbl)

####################################################################################################################################
# Cluster the signatures using Jaccard index as similarity metric
####################################################################################################################################

bp_jac_m <- jaccard(pca_loadings_tbl)

pca_ccp <- ConsensusClusterPlus::ConsensusClusterPlus(d = dist(bp_jac_m, "euclidean"), maxK = 10, reps = 1000, innerLinkage = "average", finalLinkage = "average", plot = "pdf")
pca_k <- 5
clusters <- as.factor(pca_ccp[[pca_k]]$consensusClass)

hc <- pca_ccp[[pca_k]]$consensusTree

####################################################################################################################################
# Generate consensus signatures from the 5 clusters
####################################################################################################################################

pca_gene_sigs_ranked <- lapply(levels(clusters), function(cl) {
  n <- table(clusters)[cl]
  min_n <- round(max(.25*n, 3))
  print(min_n)
  gs <- table(unlist(pca_loadings_tbl[, pca_ccp[[pca_k]]$consensusClass == cl])) %>% sort()
  gs <- gs[gs >= min_n]
  return(gs)
})
names(pca_gene_sigs_ranked) <- paste0("C", levels(clusters))

pca_gene_sigs <- lapply(pca_gene_sigs_ranked, names)

####################################################################################################################################
# Create a Gene x PB matrix containing all consensus genes to facilitate scoring
####################################################################################################################################

pca_gene <- unlist(pca_gene_sigs, recursive = F)  %>% unique() %>% unname()

pbm_list <- lapply(names(state_pb_profiles), function(st) {
  
  print(st)
  
  pbm <- state_pb_profiles[[st]][pca_gene, ]
  dim(pbm)
  
  pbm <- rowcenter(pbm)
  
  return(pbm)
})
pbm_list <- do.call(cbind, pbm_list)

####################################################################################################################################
# Score the PB profiles for the consensus signatures and create the feature_scores table
####################################################################################################################################

feature_scores <- tibble(ID = colnames(pbm_list),
                         C1 = colMeans(pbm_list[pca_gene_sigs$C1, ]),
                         C2 = colMeans(pbm_list[pca_gene_sigs$C2, ]),
                         C3 = colMeans(pbm_list[pca_gene_sigs$C3, ]),
                         C4 = colMeans(pbm_list[pca_gene_sigs$C4, ]),
                         C5 = colMeans(pbm_list[pca_gene_sigs$C5, ]),
                         Diff = C3 - C1)
feature_scores$Sample <- sapply(strsplit(feature_scores$ID, split = "_"), function(x) x[[1]][1])
feature_scores$State <- sapply(strsplit(feature_scores$ID, split = "_"), function(x) x[[2]][1])

feature_scores <- feature_scores %>%
  left_join(sample_data %>%
              select(Sample, Patient, Timepoint),
            by = "Sample")

feature_scores$Complexity <- complexity[feature_scores$ID]
feature_scores$Ncell <- n_cells[feature_scores$ID]

feature_scores$SID <- paste0(feature_scores$Patient, feature_scores$Timepoint)

####################################################################################################################################
# 3-axis coordinates for C1/3/5 scores to depict inter-tumor heterogeneity
####################################################################################################################################

feature_scores$Lineage <- apply(feature_scores[, c("C1", "C3")], 1, max)
feature_scores$Stemness <- apply(feature_scores[, c("C5", "Lineage")], 1, function(x) x[1] - x[2])

feature_scores$LineagePlot <- apply(feature_scores[, c("C1", "C3")], 1, function(x) {
  res <- max(x[1], x[2])
  if (res < 0)
    res <- runif(1, min = 0, max = 0.25)
  else if(x[1] > x[2])
    res <- -1 * res
  return (res)
})

####################################################################################################################################
# Generate feature scores per sample the and 3-axis coordinates
####################################################################################################################################

feature_scores_per_sample <- feature_scores %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(C1 = mean(C1), C2 = mean(C2), C3 = mean(C3), C4 = mean(C4), C5 = mean(C5)) %>%
  ungroup()

feature_scores_per_sample$Lineage <- apply(feature_scores_per_sample[, c("C1", "C3")], 1, max)
feature_scores_per_sample$Stemness <- apply(feature_scores_per_sample[, c("C5", "Lineage")], 1, function(x) x[1] - x[2])

feature_scores_per_sample$LineagePlot <- apply(feature_scores_per_sample[, c("C1", "C3")], 1, function(x) {
  res <- max(x[1], x[2])
  if (res < 0)
    res <- runif(1, min = 0, max = 0.25)
  else if(x[1] > x[2])
    res <- -1 * res
  return (res)
})

####################################################################################################################################
# Generate a color code for each axis to facilitate using 3 color scales
####################################################################################################################################

feature_scores_per_sample$MaxAxis <- apply(feature_scores_per_sample[, c("C1", "C3", "C5")], 1,  function(x) c("C1", "C3", "C5")[which.max(x)])

feature_scores_per_sample$C1_rescaled <- rescale(feature_scores_per_sample$C1, c(0, .3))
feature_scores_per_sample$C3_rescaled <- rescale(feature_scores_per_sample$C3, c(.35, .65))
feature_scores_per_sample$C5_rescaled <- rescale(feature_scores_per_sample$C5, c(.7, 1))

feature_scores_per_sample <- feature_scores_per_sample %>%
  mutate(axisColorCode = case_when(MaxAxis == "C1" ~ C1_rescaled,
                                   MaxAxis == "C3" ~ C3_rescaled,
                                   MaxAxis == "C5" ~ C5_rescaled,
                                   TRUE ~ NA))

axis_color_vec <- c(RColorBrewer::brewer.pal(n = 9, name = "Greens"),
                    RColorBrewer::brewer.pal(n = 9, name = "Purples"),
                    RColorBrewer::brewer.pal(n = 9, name = "Reds"))

feature_scores_per_sample$ID <- paste0(feature_scores_per_sample$Patient, feature_scores_per_sample$Timepoint)

####################################################################################################################################
# Score the cells for the 3 consensus signatures
####################################################################################################################################

m <- umi_data[, mdata$CellID]

m <- umi2upm(m)

rm <- log2(rowMeans(m) + 1)
sg <- names(rm[rm > 4])

m <- m[sg, ]

scp_scores <- lapply(unique(mdata$State), function(st) {
  
  print(st)
  
  d_st <- mdata %>%
    filter(State == st)
  
  print(nrow(d_st))
  
  m_st <- m[, d_st$CellID]
  m_st <- log2(m_st / 10 + 1)
  m_st <- as.matrix(m_st)
  
  sigScores(m = m_st, sigs = pca_gene_sigs) %>%
    as_tibble(rownames = "CellID")
})
scp_scores <- do.call(rbind, scp_scores)

scp_scores$Lineage <- apply(scp_scores[, c("C1", "C3")], 1, max)
scp_scores$Stemness <- apply(scp_scores[, c("C5", "Lineage")], 1, function(x) x[1] - x[2])

scp_scores$LineagePlot <- apply(scp_scores[, c("C1", "C3")], 1, function(x) {
  res <- max(x[1], x[2])
  if (res < 0)
    res <- runif(1, min = 0, max = 0.5)
  else if(x[1] > x[2])
    res <- -1 * res
  return (res)
})

scp_scores <- scp_scores %>%
  left_join(mdata %>%
              select(State, CellID), by = "CellID")

####################################################################################################################################
# Color-code the cells to facilitate plotting the butterflies
####################################################################################################################################

scp_scores$MaxAxis <- apply(scp_scores[, c("C1", "C3", "C5")], 1,  function(x) c("C1", "C3", "C5")[which.max(x)])

scp_scores$C1_rescaled <- rescale(scp_scores$C1, c(0, .3))
scp_scores$C3_rescaled <- rescale(scp_scores$C3, c(.35, .65))
scp_scores$C5_rescaled <- rescale(scp_scores$C5, c(.7, 1))

scp_scores <- scp_scores %>%
  mutate(axisColorCode = case_when(MaxAxis == "C1" ~ C1_rescaled,
                                   MaxAxis == "C3" ~ C3_rescaled,
                                   MaxAxis == "C5" ~ C5_rescaled,
                                   TRUE ~ NA))

####################################################################################################################################
# Generate the SCP scores for all samples (include also samples with few malignant cells)
####################################################################################################################################

d <- mdata

# Generate profiles from the 7 principal cellular states
d_state_stats <- d %>%
  filter(State %in% c("AC", "MES", "Hypoxia", "GPC", "OPC", "NPC", "Neuron")) %>%
  group_by(Patient, Timepoint, Sample, State) %>%
  summarise(n = n())

states <- c("AC", "MES", "Hypoxia", "GPC", "OPC", "NPC", "Neuron")

d_state_stats <- d_state_stats %>%
  filter(n >= 5)

# Remove "junk genes" (i.e. pseudo-genes, antisense genes etc)
junk_genes <- c(rownames(umi_data_all[[1]])[grep("\\.", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("-AS*", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("LINC", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("^RP[S|L]", rownames(umi_data_all[[1]]))],
                rownames(umi_data_all[[1]])[grep("^MT-", rownames(umi_data_all[[1]]))])

valid_genes <- rownames(umi_data_all[[1]])[rownames(umi_data_all[[1]]) %ni% junk_genes]

# Generate the pseudo-bulk profiles
state_pb_profiles_all <- lapply(unique(d_state_stats$State), function(st) {
  
  print(st)
  
  d_state_stats_tmp <- d_state_stats %>%
    filter(State == st)
  
  res <- lapply(d_state_stats_tmp$Sample, function(sname) {
    
    m <- umi_data_all[[sname]]
    
    d_tmp <- d %>%
      filter(Sample == sname, State == st)
    
    m <- m[valid_genes, d_tmp$CellID]
    
    m <- umi2upm(m)
    
    m <- log2(rowMeans(m) + 1)
    
    return(m)
  })
  res <- do.call(cbind, res)
  colnames(res) <- paste0(d_state_stats_tmp$Sample, "_", st)
  
  return(res)
})
names(state_pb_profiles_all) <- unique(d_state_stats$State)

pbm_all_list <- lapply(names(state_pb_profiles_all), function(st) {
  
  print(st)
  
  pbm <- state_pb_profiles_all[[st]][pca_gene, ]
  dim(pbm)
  
  pbm <- rowcenter(pbm)
  
  return(pbm)
})
pbm_all_list <- do.call(cbind, pbm_all_list)

feature_scores_all <- tibble(ID = colnames(pbm_all_list),
                             C1 = colMeans(pbm_all_list[pca_gene_sigs$C1, ]),
                             C3 = colMeans(pbm_all_list[pca_gene_sigs$C3, ]),
                             C5 = colMeans(pbm_all_list[pca_gene_sigs$C5, ]))
feature_scores_all$Sample <- sapply(strsplit(feature_scores_all$ID, split = "_"), function(x) x[[1]][1])
feature_scores_all$State <- sapply(strsplit(feature_scores_all$ID, split = "_"), function(x) x[[2]][1])

feature_scores_all <- feature_scores_all %>%
  left_join(sample_data %>%
              select(Sample, Patient, Timepoint),
            by = "Sample")

feature_scores_all_per_sample <- feature_scores_all %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(C1 = mean(C1), C3 = mean(C3), C5 = mean(C5)) %>%
  ungroup()

feature_scores_all_per_sample$Lineage <- apply(feature_scores_all_per_sample[, c("C1", "C3")], 1, max)
feature_scores_all_per_sample$Stemness <- apply(feature_scores_all_per_sample[, c("C5", "Lineage")], 1, function(x) x[1] - x[2])

feature_scores_all_per_sample$LineagePlot <- apply(feature_scores_all_per_sample[, c("C1", "C3")], 1, function(x) {
  res <- max(x[1], x[2])
  if (res < 0)
    res <- runif(1, min = 0, max = 0.25)
  else if(x[1] > x[2])
    res <- -1 * res
  return (res)
})

feature_scores_all_per_sample <- feature_scores_all_per_sample %>%
  group_by(Sample) %>%
  mutate(SCP = c("SCP-ECM", "SCP-Neuronal", "SCP-Glial")[which.max(c(C1, C3, C5))]) %>%
  ungroup()
