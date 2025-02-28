#################################
## Title: Extended Data Figure 4 panels e-f in Nomura et al - contribution of lab and medical center per meta-program
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a barplot showing the contribution of data-generating lab and medical center per meta-program
#################################

####################################################################################################################################
# Figure EDF4e - Lab contribution per MP
####################################################################################################################################

nmf_clusters <- malignant_nmf_metaprograms$clusters

res <- lapply(names(nmf_clusters), function(xname) {
  x <- nmf_clusters[[xname]]
  res <- sapply(strsplit(unique(x), split = "_"), function(y) y[[1]][1])
  table(sample_data$Lab[sample_data$Sample %in% res]) %>%
    melt() %>%
    mutate(cluster = xname)
})
res <- do.call(rbind, res) %>%
  as_tibble()

res$Var1 <- as.character(res$Var1)

n_samples_per_lab <- meta_data %>%
  select(Sample, Lab) %>%
  filter(!duplicated(Sample)) %>%
  group_by(Lab) %>%
  summarise(n = n(), .groups = "drop")

n_samples_per_lab <- setNames(n_samples_per_lab$n, n_samples_per_lab$Lab)

res$value <- res$value / n_samples_per_lab[res$Var1]

res$MP <- gsub("Cluster", "MP", res$cluster)

res <- res %>%
  filter(MP != "MP_16")

res$MP <- factor(res$MP, paste0("MP_", 1:15))

ggplot(res, aes(x = MP, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_discrete(name = "Lab") +
  scale_y_continuous(labels = percent) +
  xlab("Meta-program") +
  ylab("% contributing samples\n(proportional to lab-specific samples)") +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 12)) +
  theme(panel.grid.major = element_line())

####################################################################################################################################
# Figure EDF4f - MC contribution per MP
####################################################################################################################################

nmf_clusters <- malignant_nmf_metaprograms$clusters

res <- lapply(names(nmf_clusters), function(xname) {
  x <- nmf_clusters[[xname]]
  res <- sapply(strsplit(unique(x), split = "_"), function(y) y[[1]][1])
  table(sample_data$Origin[sample_data$Sample %in% res]) %>%
    melt() %>%
    mutate(cluster = xname)
})
res <- do.call(rbind, res) %>%
  as_tibble()

res$Var1 <- as.character(res$Var1)

n_samples_per_origin <- meta_data %>%
  select(Sample, Origin) %>%
  filter(!duplicated(Sample)) %>%
  group_by(Origin) %>%
  summarise(n = n(), .groups = "drop")

n_samples_per_origin <- setNames(n_samples_per_origin$n, n_samples_per_origin$Origin)

res$value <- res$value / n_samples_per_origin[res$Var1]

res$MP <- gsub("Cluster", "MP", res$cluster)

res <- res %>%
  filter(MP != "MP_16")

res$MP <- factor(res$MP, paste0("MP_", 1:15))

ggplot(res, aes(x = MP, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_discrete(name = "Origin") +
  scale_y_continuous(labels = percent) +
  xlab("Meta-program") +
  ylab("% contributing samples\n(proportional to institution-specific samples)") +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 12)) +
  theme(panel.grid.major = element_line())
