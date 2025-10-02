#################################
## Title: Extended Data Figure 5 panel d in Nomura et al - number of detected genes in hybrid vs. singular cells
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a boxplot comparing the number of detected genes in hybrid vs. singular cells
#################################

d <- mdata %>%
  mutate(isHybrid = mdata$HybridID != "Singular")

dm <- d %>%
  group_by(Sample, isHybrid) %>%
  summarise(comp = mean(nFeature_RNA), SD = sd(nFeature_RNA), N = n(), SE = SD / sqrt(N), .groups = "drop")

ggboxplot(data = dm, x = "isHybrid", y = "comp", add = "jitter") +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

t.test(x = dm$comp[d$isHybrid], y = dm$comp[!d$isHybrid])
