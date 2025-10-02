#################################
## Title: Figure 5 panel d in Nomura et al - multi-layer graph
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Generate a graph showing relationship between transriptomic and genetic layers
#################################

####################################################################################################################################
# Basic pairs (Comp-Mal, Comp-BP, Mal-BP)
####################################################################################################################################

feature_scores_all_per_sample$Diff <- apply(feature_scores_all_per_sample[, c("C1", "C3", "C5")] %>% as.matrix(), 1, function(x) sort(x)[3] - sort(x)[2])
feature_scores_all_per_sample$MaxScore <- apply(feature_scores_all_per_sample[, c("C1", "C3", "C5")] %>% as.matrix(), 1, max)

feature_scores_all_per_sample$SCP[feature_scores_all_per_sample$MaxScore < 0 | feature_scores_all_per_sample$Diff < .25] <- "Mixed"

feature_scores_all_per_sample$SCP <- factor(feature_scores_all_per_sample$SCP, c("SCP-ECM", "SCP-Neuronal", "SCP-Glial", "Mixed"))

gbm_subtypes_tbl$SCP <- setNames(feature_scores_all_per_sample$SCP, feature_scores_all_per_sample$Sample)[gbm_subtypes_tbl$Sample]

d <- gbm_subtypes_tbl %>%
  filter(!is.na(MalCluster))

d$Pair1 <- paste0(d$CompCluster, "_", d$MalCluster)
d$Pair2 <- paste0(d$CompCluster, "_", d$SCP)
d$Pair3 <- paste0(d$MalCluster, "_", d$SCP)

pair1_obs_vs_exp <- lapply(unique(d$Pair1), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  nobs <- sum(d$Pair1 == p)
  prob <- (sum(d$CompCluster == p1) / nrow(d)) * (sum(d$MalCluster == p2) / nrow(d))
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / nrow(d), fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = nrow(d), p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair1_obs_vs_exp <- do.call(rbind, pair1_obs_vs_exp)

pair2_obs_vs_exp <- lapply(unique(d$Pair2), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  nobs <- sum(d$Pair2 == p)
  prob <- (sum(d$CompCluster == p1) / nrow(d)) * (sum(d$SCP == p2) / nrow(d))
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / nrow(d), fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = nrow(d), p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair2_obs_vs_exp <- do.call(rbind, pair2_obs_vs_exp)

pair3_obs_vs_exp <- lapply(unique(d$Pair3), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  nobs <- sum(d$Pair3 == p)
  prob <- (sum(d$MalCluster == p1) / nrow(d)) * (sum(d$SCP == p2) / nrow(d))
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / nrow(d), fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = nrow(d), p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair3_obs_vs_exp <- do.call(rbind, pair3_obs_vs_exp)

sig_pairs_tbl <- rbind(pair1_obs_vs_exp %>%
                         filter(pval < .1),
                       pair2_obs_vs_exp %>%
                         filter(pval < .05),
                       pair3_obs_vs_exp %>%
                         filter(pval < .05))

####################################################################################################################################
# Add CNA drivers
####################################################################################################################################

driver_data <- readRDS(paste0(DATA_ROOT, "driver_cna_status_care_synapse_20231109.RDS")) %>%
  select(-aliquot_barcode, -Patient, -Timepoint)

drivers <- colnames(driver_data)[2:6]

driver_data <- driver_data %>%
  left_join(gbm_subtypes_tbl %>%
              select(ID, CompCluster, MalCluster, SCP), by = "ID")

d <- driver_data %>%
  filter(!is.na(MalCluster), !is.na(SCP))

N <- length(unique(d$ID))

d <- melt(d, variable.name = "Gene", value.name = "CNA") %>%
  as_tibble() %>%
  filter(CNA == 1) %>%
  mutate(Gene = as.character(Gene))

d$Pair1 <- paste0(d$Gene, "_", d$CompCluster)
d$Pair2 <- paste0(d$Gene, "_", d$MalCluster)
d$Pair3 <- paste0(d$Gene, "_", d$SCP)

pair1_obs_vs_exp <- lapply(unique(d$Pair1), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% driver_data$ID, CompCluster == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair1 == p & d$CNA == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair1_obs_vs_exp <- do.call(rbind, pair1_obs_vs_exp)

pair2_obs_vs_exp <- lapply(unique(d$Pair2), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% driver_data$ID, MalCluster == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair2 == p & d$CNA == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair2_obs_vs_exp <- do.call(rbind, pair2_obs_vs_exp)

pair3_obs_vs_exp <- lapply(unique(d$Pair3), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% driver_data$ID, SCP == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair3 == p & d$CNA == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair3_obs_vs_exp <- do.call(rbind, pair3_obs_vs_exp)

sig_pairs_tbl <- rbind(sig_pairs_tbl,
                       rbind(pair1_obs_vs_exp %>%
                               filter(pval < .1),
                             pair2_obs_vs_exp %>%
                               filter(pval < .1),
                             pair3_obs_vs_exp %>%
                               filter(pval < .1)))

####################################################################################################################################
# Add SNVs
####################################################################################################################################

snv_data <- readRDS(paste0(DATA_ROOT, "driver_snv_status_care_synapse_20231109.RDS")) %>%
  select(-aliquot_barcode, -case_barcode)

snvs <- colnames(snv_data)[2:10]

snv_data <- snv_data %>%
  left_join(gbm_subtypes_tbl %>%
              select(ID, CompCluster, MalCluster, SCP), by = "ID")

d <- snv_data %>%
  filter(!is.na(MalCluster), !is.na(SCP))
dim(d)

N <- length(unique(d$ID))

d <- melt(d, variable.name = "Gene", value.name = "SNV") %>%
  as_tibble() %>%
  filter(SNV == 1) %>%
  mutate(Gene = as.character(Gene))

d$Pair1 <- paste0(d$Gene, "_", d$CompCluster)
d$Pair2 <- paste0(d$Gene, "_", d$MalCluster)
d$Pair3 <- paste0(d$Gene, "_", d$SCP)

pair1_obs_vs_exp <- lapply(unique(d$Pair1), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% snv_data$ID, CompCluster == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair1 == p & d$SNV == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair1_obs_vs_exp <- do.call(rbind, pair1_obs_vs_exp)

pair2_obs_vs_exp <- lapply(unique(d$Pair2), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% snv_data$ID, MalCluster == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair2 == p & d$SNV == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair2_obs_vs_exp <- do.call(rbind, pair2_obs_vs_exp)

pair3_obs_vs_exp <- lapply(unique(d$Pair3), function(p) {
  p1 <- strsplit(x = p, split = "_")[[1]][1]
  p2 <- strsplit(x = p, split = "_")[[1]][2]
  
  n2 <- gbm_subtypes_tbl %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% snv_data$ID, SCP == p2) %>%
    nrow()
  
  nobs <- sum(d$Pair3 == p & d$SNV == 1)
  prob <- (sum(d$Gene == p1) / N) * (n2 / N)
  
  res <- tibble(Pair = p, nobs = nobs, fobs = nobs / N, fexp = prob, FC = fobs / fexp, log2FC = log2(FC),
                pval = binom.test(x = nobs, n = N, p = prob, alternative = "greater")$p.value)
  
  return(res)
})
pair3_obs_vs_exp <- do.call(rbind, pair3_obs_vs_exp)

sig_pairs_tbl <- rbind(sig_pairs_tbl,
                       rbind(pair1_obs_vs_exp %>%
                               filter(pval < .05),
                             pair2_obs_vs_exp %>%
                               filter(pval < .05),
                             pair3_obs_vs_exp %>%
                               filter(pval < .05)))

####################################################################################################################################
# Generate graph and plot
####################################################################################################################################

edges <- lapply(sig_pairs_tbl$Pair, function(h) {
  
  h <- strsplit(x = as.character(h), split = "_")
  
  vs1 <- h[[1]][1]
  vs2 <- h[[1]][2]
  data.frame(S1 = vs1, S2 = vs2)
})
edges <- do.call(rbind, edges)
edges$log2FC <- sig_pairs_tbl$log2FC
edges$FC <- sig_pairs_tbl$FC

pairs_graph <- graph_from_data_frame(edges, directed = F)

colors <- RColorBrewer::brewer.pal(6, "Spectral")

E(pairs_graph)$weight <- edges$FC
E(pairs_graph)$arrow.size <- edges$FC
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% levels(gbm_subtypes_tbl$SCP)] <- colors[1]
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% levels(gbm_subtypes_tbl$MalCluster)] <- colors[2]
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% levels(gbm_subtypes_tbl$CompCluster)] <- colors[3]
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% drivers] <- colors[4]
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% snvs] <- colors[5]
V(pairs_graph)$vertex.color[rownames(pairs_graph[]) %in% cnas] <- colors[6]

plot(pairs_graph,
     edge.width = E(pairs_graph)$weight,
     vertex.color = V(pairs_graph)$vertex.color, vertex.label.color="black",
     edge.label.cex = 1,
     edge.curved = F,
     vertex.size = 10, vertex.label.cex = 1.5)

multi_layer_groups <- edges

####################################################################################################################################
# Classify samples to ecosystems
####################################################################################################################################

ecosystems_vec <- components(pairs_graph)$membership
ecosystems_vec <- setNames(paste0("ES", ecosystems_vec), names(ecosystems_vec))

d <- gbm_subtypes_tbl %>%
  select(Sample, ID, CompCluster, MalCluster, SCP) %>%
  mutate(ID = as.character(ID),
         CompCluster = as.character(CompCluster),
         MalCluster = as.character(MalCluster),
         SCP = as.character(SCP))

d <- d %>%
  filter(!is.na(MalCluster), !is.na(SCP))

d <- d %>%
  melt(id.vars = c("Sample", "ID")) %>%
  as_tibble() %>%
  mutate(variable = as.character(variable))

d$ES_comp <- ecosystems_vec[d$value]

d <- d %>%
  filter(!is.na(ES_comp))

d <- d %>%
  group_by(Sample, ID, ES_comp) %>%
  mutate(ID = as.character(ID)) %>%
  summarise(n = n(), .groups = "drop")

es_class <- lapply(unique(d$ID), function(id) {
  x <- d %>%
    filter(ID == id) %>%
    arrange(desc(n))
  
  es <- NA
  
  if(nrow(x) == 1 & max(x$n) > 1) {
    es <- x$ES_comp
  } else if(nrow(x) > 1) {
    if(x$n[1] > 2 & x$n[2] <= 1) {
      es <- x$ES_comp[1]
    } else if((x$n[1] == 2 & x$n[2] == 0) | (x$ES_comp[1] == "ES2" & x$n[1] == 2))  {
      es <- x$ES_comp[1]
    }
  }
  
  return(tibble(Sample = x$Sample[1], ID = x$ID[1], ES = es, n = max(x$n)))
})
es_class <- do.call(rbind, es_class)

es_vec <- setNames(rep(NA, nrow(gbm_subtypes_tbl)), gbm_subtypes_tbl$Sample)

es_vec[es_class$Sample] <- es_class$ES

gbm_subtypes_tbl$ES <- es_vec[gbm_subtypes_tbl$Sample]

d <- gbm_subtypes_tbl %>%
  mutate(ES = case_when(CompCluster == "Mixed" | MalCluster == "Mixed" | SCP == "Mixed" ~ NA,
                        TRUE ~ ES))

####################################################################################################################################
# Statistical test for ecosystems
####################################################################################################################################

d <- gbm_subtypes_tbl %>%
  select(Sample, ID, CompCluster, MalCluster, SCP) %>%
  mutate(ID = as.character(ID),
         CompCluster = as.character(CompCluster),
         MalCluster = as.character(MalCluster),
         SCP = as.character(SCP))

d <- d %>%
  filter(!is.na(MalCluster), !is.na(SCP))

d <- d %>%
  melt(id.vars = c("Sample", "ID")) %>%
  as_tibble() %>%
  mutate(variable = as.character(variable))

perm_test <- lapply(1:25, function(i) {
  print(i)
  
  d_shuff <- d %>%
    group_by(variable) %>%
    mutate(value = sample(x = value, size = n())) %>%
    ungroup()
  
  d_shuff$ES_comp <- ecosystems_vec[d_shuff$value]
  
  d_shuff <- d_shuff %>%
    filter(!is.na(ES_comp))
  
  d_shuff <- d_shuff %>%
    group_by(Sample, ID, ES_comp) %>%
    mutate(ID = as.character(ID)) %>%
    summarise(n = n(), .groups = "drop")
  
  es_class <- lapply(unique(d_shuff$ID), function(id) {
    x <- d_shuff %>%
      filter(ID == id) %>%
      arrange(desc(n))
    
    es <- NA
    
    if(nrow(x) == 1 & max(x$n) > 1) {
      es <- x$ES_comp
    } else if(nrow(x) > 1) {
      if(x$n[1] > 2 & x$n[2] <= 1) {
        es <- x$ES_comp[1]
      } else if((x$n[1] == 2 & x$n[2] == 0) | (x$ES_comp[1] == "ES2" & x$n[1] == 2))  {
        es <- x$ES_comp[1]
      }
    }
    
    return(tibble(Sample = x$Sample[1], ID = x$ID[1], ES = es, n = max(x$n)))
  })
  es_class <- do.call(rbind, es_class)
  
  es_vec <- setNames(rep(NA, nrow(gbm_subtypes_tbl)), gbm_subtypes_tbl$Sample)
  
  es_vec[es_class$Sample] <- es_class$ES
  as_tibble(table(es_vec))
})
perm_test <- do.call(rbind, perm_test)

perm_test_sum <- perm_test %>%
  group_by(es_vec) %>%
  summarise(n = mean(n))

fisher.test(matrix(data = c(perm_test_sum$n, table(gbm_subtypes_tbl$ES)), nrow = 2, byrow = T))

fisher.test(matrix(data = c(c(sum(perm_test_sum$n), 120 - sum(perm_test_sum$n)),
                            c(sum(table(gbm_subtypes_tbl$ES)), 120 - sum(table(gbm_subtypes_tbl$ES)))), nrow = 2, byrow = T))

N <- gbm_subtypes_tbl %>%
  filter(!is.na(MalCluster), !is.na(SCP)) %>%
  nrow()
nobs <- table(gbm_subtypes_tbl$ES)
nexp <- round(perm_test_sum$n, 0)

sapply(1:3, function(i) binom.test(x = nobs[i], n = N, p = nexp[i] / N)$p.value)
