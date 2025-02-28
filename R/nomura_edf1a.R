#################################
## Title: Extended Data Figure 1 panel a in Nomura et al - quality control metrics
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Quality control metrics - total counts per cell, number of detected genes per cell and mitochondrial expression
#################################

d <- meta_data %>%
  as_tibble()

ggplot(d, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100, color = "black", fill = "dodgerblue", alpha = .25) +
  geom_vline(xintercept = quantile(d$nCount_RNA, .1), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = quantile(d$nCount_RNA, .9), linetype = "dashed", color = "black", size = 1) +
  xlab("Total counts per cell [log10]") +
  ylab("Number of cells") +
  scale_x_log10() +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

ggplot(d, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, color = "black", fill = "dodgerblue", alpha = .25) +
  geom_vline(xintercept = quantile(d$nFeature_RNA, .1), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = quantile(d$nFeature_RNA, .9), linetype = "dashed", color = "black", size = 1) +
  xlab("Number of detected genes per cell") +
  ylab("Number of cells") +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

ggplot(d, aes(x = Percent.MT)) +
  geom_histogram(bins = 100, color = "black", fill = "dodgerblue", alpha = .25) +
  xlab("Mitochondrial expression [%]") +
  ylab("Number of cells") +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())
