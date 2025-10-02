#################################
## Title: Extended Data Figure 2d in Spitzer et al - baseline porifle differences across timepoints
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a boxplot of baseline differences across timepoints
#################################

feature_scores_all_per_sample <- feature_scores_all_per_sample %>%
  left_join(gbm_subtypes_tbl %>%
              select(Sample, PvsR), by = "Sample")

d <- feature_scores_all_per_sample %>%
  dplyr::select(Sample, PvsR, C1, C3, C5) %>%
  melt() %>%
  as_tibble() %>%
  dplyr::rename(BP = variable, Score = value) %>%
  mutate(BP = as.character(BP),
         BP = case_when(BP == "C1" ~ "SCP-ECM",
                        BP == "C3" ~ "SCP-Neuronal",
                        BP == "C5" ~ "SCP-Glial"),
         BP = factor(BP, c("SCP-ECM", "SCP-Neuronal", "SCP-Glial")))

ggboxplot(data = d, x = "BP", y = "Score", color = "PvsR", add = "jitter") +
  scale_color_brewer(name = "", palette = "Set2") +
  xlab("") +
  ylab("Mean score") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "right", axis.text.x = element_text(size = 12, angle = 90)) +
  theme(panel.grid.major = element_line(), strip.text = element_text(size = 16))

d <- rbind(d,
           d %>%
             group_by(Sample, PvsR) %>%
             summarise(Score = Score[BP == "SCP-Neuronal"] - Score[BP == "SCP-Glial"],
                       BP = "NGdiff", .groups = "drop") %>%
             dplyr::select(Sample, PvsR, BP, Score))

d_stats <- d %>%
  group_by(BP, PvsR) %>%
  summarise(Mean = mean(Score), SD = sd(Score), n = n(), SE = SD / sqrt(n),
            .groups = "drop")

pval_tbl <- d %>%
  group_by(BP) %>%
  summarise(pval = t.test(Score[PvsR == "Primary"], Score[PvsR != "Primary"], var.equal = T)$p.value,
            .groups = "drop")
