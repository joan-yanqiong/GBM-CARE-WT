#################################
## Title: Extended Data Figure 5a in Spitzer et al - interval between first and second resection stratified by MGMT methylation status
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce survival curves for the surgical interval stratified by MGMT methylation status
#################################

d <- mgmt_md %>%
  mutate(one = 1) %>%
  filter(MGMT != "NOS", Timepoint == "T1") %>%
  mutate(MGMT = factor(MGMT, c("MET", "UM", "NOS")))

km_fit <- survfit(Surv(SurgicalInterval, one) ~ MGMT, data = d)
ggsurvplot(fit = km_fit, data = d,
           pval = T,
           risk.table = T,
           # conf.int = T,
           # surv.median.line = "hv",
           break.x.by = 10,    
           xlab = "Time in months",
           xlim = c(0, 30),
           # fontsize = 10,
           size = 1,
           # risk.table.fontsize = 10,
           legend.title = "MGMT status",
           # legend.labs = c("Other", "Glio-neural"),
           # palette = "Set1",
           palette = c("red", "black", "grey"),
           surv.scale = "percent", ggtheme = theme_bw()) #+
