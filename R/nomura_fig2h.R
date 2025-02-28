#################################
## Title: Figure 2 panel h in Nomura et al - butterfly plot emphasizing the GPC state
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a butterfly plot (similar to Neftel et al) emphasizing the GPC state
#################################

# include all samples with at least 100 cells in both timepoints
valid_patients <- mdata %>%
  group_by(Patient, Timepoint) %>%
  summarise(n = n()) %>%
  arrange(n) %>%
  filter(n >= 100) %>%
  filter("T1" %in% Timepoint & "T2" %in% Timepoint) %>%
  pull(Patient) %>%
  unique()

# Donw-sample each sample to 100 cells for easier plotting
d <- mdata %>%
  filter(Patient %in% valid_patients) %>%
  group_by(Patient, Timepoint) %>%
  sample_n(100)

knn_res <- FNN::get.knn(d[, c("Dx", "Dy")], k = 10)
d$isGPC <- d$State == "GPC"
d$GPC_density <- sapply(1:nrow(knn_res$nn.index), function(i) sum(d$isGPC[knn_res$nn.index[i, ]]))
table(d$GPC_density)

d$GPC_density <- d$GPC_density / 10

d <- d %>% arrange(GPC_density)

ggplot(data = d, mapping = aes(x = Dx, y = Dy)) +
  geom_point(mapping = aes(fill = GPC_density), shape = 21, color = "black", size = 4,
             show.legend = c("color" = TRUE, "fill" = TRUE, "alpha" = FALSE)) +
  scale_fill_gradient2(name = "% GPC neighbors", low = "black", mid = "brown", high = "red", oob = squish, limits = c(0, .4), midpoint = .2, labels = c("0%", "", "20%", "", "40%")) +
  xlab("Relative meta-module score\n[log2(|SC1 - SC2| + 1)]") +
  ylab("Relative meta-module score\n[log2(|SC1 - SC2| + 1)]") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr() +
  theme(panel.grid.major = element_line())
