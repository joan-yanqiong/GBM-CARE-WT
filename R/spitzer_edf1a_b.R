#################################
## Title: Extended Data Figure 1a-b in Spitzer et al - number of cells and average number of detected genes across timepoints
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce boxplots of number of cells and average number of detected genes across timepoints
#################################

####################################################################################################################################
# Number of cells per sample across timepoints
####################################################################################################################################

meta_data %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(Ncells = n()) %>%
  mutate(Timepoint = factor(Timepoint, c("T1", "T2", "T3"))) %>%
  ggboxplot(x = "Timepoint", y = "Ncells", add = "jitter", fill = "Timepoint") +
  scale_color_brewer(name = "", palette = "Set1") +
  xlab("Timepoint") +
  ylab("Number of cells\nper sample") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "none") +
  theme(panel.grid.major = element_line(), strip.text = element_text(size = 16))

meta_data %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(Ncells = n()) %>%
  ungroup() %>%
  summarise(p1 = wilcox.test(x = Ncells[Timepoint == "T1"], y = Ncells[Timepoint == "T2"])$p.value,
            p2 = wilcox.test(x = Ncells[Timepoint == "T1"], y = Ncells[Timepoint == "T3"])$p.value,
            p3 = wilcox.test(x = Ncells[Timepoint == "T2"], y = Ncells[Timepoint == "T3"])$p.value)

####################################################################################################################################
# Average number of detected genes per sample across timepoints
####################################################################################################################################

meta_data %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(Comp = mean(nFeature_RNA)) %>%
  mutate(Timepoint = factor(Timepoint, c("T1", "T2", "T3"))) %>%
  ggboxplot(x = "Timepoint", y = "Comp", add = "jitter", fill = "Timepoint") +
  scale_color_brewer(name = "", palette = "Set1") +
  xlab("Timepoint") +
  ylab("Average number of detected\ngenes per sample") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "none") +
  theme(panel.grid.major = element_line(), strip.text = element_text(size = 16))

meta_data %>%
  group_by(Sample, Patient, Timepoint) %>%
  summarise(Comp = mean(nFeature_RNA)) %>%
  ungroup() %>%
  summarise(p1 = wilcox.test(x = Comp[Timepoint == "T1"], y = Comp[Timepoint == "T2"])$p.value,
            p2 = wilcox.test(x = Comp[Timepoint == "T1"], y = Comp[Timepoint == "T3"])$p.value,
            p3 = wilcox.test(x = Comp[Timepoint == "T2"], y = Comp[Timepoint == "T3"])$p.value)
