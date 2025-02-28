#################################
## Title: Figure 2 panel g in Nomura et al - lineage model inferred from hybrdis
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a graph showing the connectivity between states
#################################

####################################################################################################################################
# Compute the observed vs. expected proportion of hybrids
####################################################################################################################################

obs_hybrid_freq <- mdata %>%
  group_by(HybridID) %>%
  summarise(n = n()) %>%
  mutate(N = sum(n), Freq = n / N)

exp_hybrid_freq <- mdata %>%
  group_by(State) %>%
  summarise(n = n()) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  arrange(State)

hybrid_ids <- unique(mdata$HybridID)

exp_hybrid <- lapply(hybrid_ids, function(hid) {
  hid1 <- strsplit(hid, split = "-")[[1]][1]
  hid2 <- strsplit(hid, split = "-")[[1]][2]
  res <- tibble(Hybrid = hid,
                ExpFreq = exp_hybrid_freq$Freq[exp_hybrid_freq$State == hid1] * exp_hybrid_freq$Freq[exp_hybrid_freq$State == hid2])
  return(res)
})
exp_hybrid <- do.call(rbind, exp_hybrid)

obs_hybrids <- melt(table(mdata$HybridID) / nrow(mdata)) %>%
  as_tibble()
colnames(obs_hybrids) <- c("Hybrid", "ObsFreq")

obs_vs_exp_hybrids <- exp_hybrid %>%
  left_join(obs_hybrids, by = "Hybrid")

obs_vs_exp_hybrids <- obs_vs_exp_hybrids[!is.na(obs_vs_exp_hybrids$ExpFreq) & !is.na(obs_vs_exp_hybrids$ObsFreq), ]

obs_vs_exp_hybrids$FC <- obs_vs_exp_hybrids$ObsFreq / obs_vs_exp_hybrids$ExpFreq
obs_vs_exp_hybrids$log2FC <- log2(obs_vs_exp_hybrids$FC)

obs_vs_exp_hybrids$Hybrid <- factor(obs_vs_exp_hybrids$Hybrid, obs_vs_exp_hybrids %>%
                                      arrange(log2FC) %>% pull(Hybrid))

technical_factor <- obs_vs_exp_hybrids %>%
  filter(log2FC < -1) %>%
  summarise(x = mean(FC)) %>%
  pull(x)

obs_vs_exp_hybrids$ExpFreqTechnical <- obs_vs_exp_hybrids$ExpFreq * technical_factor

obs_vs_exp_hybrids$log2FC_2 <- log2(obs_vs_exp_hybrids$ObsFreq / obs_vs_exp_hybrids$ExpFreqTechnical)

obs_vs_exp_hybrids$N <- obs_hybrid_freq$N[1]
obs_vs_exp_hybrids$NexpTechnical <- round(obs_vs_exp_hybrids$ExpFreqTechnical * obs_vs_exp_hybrids$N, 0)
obs_vs_exp_hybrids$Nexp <- round(obs_vs_exp_hybrids$ExpFreq * obs_vs_exp_hybrids$N, 0)
obs_vs_exp_hybrids$Nobs <- obs_vs_exp_hybrids$ObsFreq * obs_vs_exp_hybrids$N

obs_vs_exp_hybrids$Nobs_per1000 <- round(obs_vs_exp_hybrids$ObsFreq * 1000, 0)
obs_vs_exp_hybrids$Nexp_per1000 <- round(obs_vs_exp_hybrids$ExpFreq * 1000, 0)
obs_vs_exp_hybrids$NexpTechnical_per1000 <- round(obs_vs_exp_hybrids$ExpFreqTechnical * 1000, 0)

obs_vs_exp_hybrids <- obs_vs_exp_hybrids %>%
  group_by(Hybrid) %>%
  mutate(pval = fisher.test(x = matrix(data = c(c(Nobs_per1000, 1000 - Nobs_per1000),
                                                c(NexpTechnical_per1000, 1000 - NexpTechnical_per1000)),
                                       byrow = T, nrow = 2), alternative = "greater")$p.value) %>%
  ungroup()

obs_vs_exp_hybrids$padj <- p.adjust(obs_vs_exp_hybrids$pval, "fdr")

####################################################################################################################################
# Compute and plot the graph
####################################################################################################################################

d <- obs_vs_exp_hybrids %>%
  filter(sigRes)

edges <- lapply(d$Hybrid, function(h) {
  
  h <- strsplit(x = as.character(h), split = "-")
  
  vs1 <- h[[1]][1]
  vs2 <- h[[1]][2]
  data.frame(S1 = vs1, S2 = vs2)
})
edges <- do.call(rbind, edges)
edges$Freq <- d$ObsFreq
edges$log2FC <- d$log2FC_2
edges$FC <- 2^edges$log2FC

edges <- edges %>%
  filter(S2 != "Stress")

hybrid_graph <- graph_from_data_frame(edges, directed = F)

E(hybrid_graph)$weight <- edges$log2FC
E(hybrid_graph)$arrow.size <- edges$log2FC

plot(hybrid_graph,
     edge.arrow.size = E(hybrid_graph)$arrow.size, edge.width = E(hybrid_graph)$weight, edge.label.cex = 2,
     edge.label = E(hybrid_graph)$label, edge.curved = F,
     vertex.size = 25, vertex.label.cex = 2)

