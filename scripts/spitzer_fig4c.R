#################################
## Title: Spitzer Figure 4 panel c - interval between the two surgeries stratified by MGMT exression level
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate survival curves for the surgical interval stratified by MGMT expression level
#################################

d <- mgmt_md %>%
  filter(Timepoint == "T1") %>%
  mutate(one = 1)

km_fit <- survfit(Surv(SurgicalInterval, one) ~ MGMT_exp_level, data = d)

ggsurvplot(fit = km_fit, data = d,
           pval = T,
           risk.table = T,
           xlim = c(0, 30),
           break.x.by = 5,
           palette = c("black", "grey", "red"))
