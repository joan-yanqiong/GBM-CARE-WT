#################################
## Title: Extended Data Figure 6 panel d in Nomura et al - observed vs. expected within-sample variability in BP-Neuronal - BP-ECM difference
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: observed vs. expected within-sample variability in BP-Neuronal - BP-ECM difference
#################################

obs_sample_sd <- feature_scores %>%
  group_by(Sample) %>%
  summarise(n = n(), SD = sd(Diff)) %>%
  filter(n >= 2) %>%
  ungroup()

exp_sample_sd <- lapply(1:100, function(i) {
  res <- feature_scores %>%
    sample_n(3) %>%
    summarise(n = n(), SD = sd(Diff))
  return(res)
})
exp_sample_sd <- do.call(rbind, exp_sample_sd)

obs_vs_exp_sample_sd <- rbind(obs_sample_sd %>% mutate(DataType = paste0("Observed")) %>% select(SD, DataType),
                              exp_sample_sd %>% select(SD) %>% mutate(DataType = "Expected"))

ggplot(obs_vs_exp_sample_sd, aes(x = DataType, y = SD)) + 
  ggdist::stat_halfeye(aes(fill = DataType), color = "black", adjust = .5, width = .6, justification = -.2, .width = 0, point_colour = NA) + 
  geom_boxplot(width = .12, outlier.color = NA) +
  gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .3) +
  coord_cartesian(xlim = c(1.2, NA)) +
  scale_fill_brewer(name = "", palette = "Set1") +
  xlab("") +
  ylab("BPP-Neuronal - BP-ECM\ndifference intra-sample SD") +
  scale_y_continuous(limits = c(0, 3), oob = squish) +
  theme_gbm_pvsr() +
  theme(panel.grid.major = element_line())
