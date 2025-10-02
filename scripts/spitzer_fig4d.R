#################################
## Title: Spitzer Figure 4 panel d - interval between the two surgeries stratified by MGMT exression level in the MGMT methylated group
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate survival curves for the surgical interval stratified by MGMT expression level in the MGMT methylated group
#################################

d <- mgmt_md %>%
  filter(Timepoint == "T1", MGMT == "MET") %>%
  mutate(MGMT_exp_level_bin = ifelse(MGMT_exp_level == "Low", "Low", "Int/High")) %>%
  mutate(one = 1)

km_fit <- survfit(Surv(SurgicalInterval, one) ~ MGMT_exp_level_bin, data = d)

ggsurvplot(fit = km_fit, data = d,
           pval = T,
           risk.table = T,
           xlim = c(0, 30),
           break.x.by = 5,
           