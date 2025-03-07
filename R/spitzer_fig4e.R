#################################
## Title: Spitzer Figure 4 panel e - interval between the two surgeries stratified by MGMT methylation status in the MGMT-high group
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate survival curves for the surgical interval stratified by MGMT methylation status in the MGMT-high group
#################################

d <- mgmt_md %>%
  filter(Timepoint == "T1") %>%
  filter(MGMT_exp_level != "Low") %>%
  mutate(one = 1)

km_fit <- survfit(Surv(SurgicalInterval, one) ~ MGMT, data = d)

ggsurvplot(fit = km_fit, data = d,
           pval = T,
           risk.table = T,
           xlim = c(0, 30),
           break.x.by = 5,
           palette = c("red", "grey", "black"))
