#################################
## Title: Extended Data Figure 2g in Spitzer et al - survival curves
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce survival curves for composition and malignant state groups
#################################

####################################################################################################################################
# Composition group
####################################################################################################################################

d <- gbm_subtypes_tbl %>%
  filter(Timepoint == "T2")

km_fit <- survfit(Surv(SurgicalInterval2, VT) ~ CompCluster, data = d)
ggsurvplot(fit = km_fit, data = d,
           pval = T,
           risk.table = T,
           break.x.by = 10,
           xlab = "Time in months",
           size = 1,
           legend.title = "Composition group",
           palette = "Set1", surv.scale = "percent", ggtheme = theme_bw()) #+

####################################################################################################################################
# Malignant state group
####################################################################################################################################

d <- gbm_subtypes_tbl %>%
  filter(Timepoint == "T2") %>%
  mutate(MalCluster = as.character(MalCluster))

d$MalCluster[is.na(d$MalCluster)] <- "Unknown"

km_fit <- survfit(Surv(SurgicalInterval2, VT) ~ MalCluster, data = d)
ggsurvplot(fit = km_fit, data = d,
           pval = T,
           risk.table = T,
           xlim = c(0, 30),
           break.x.by = 10,    
           xlab = "Time in months",
           size = 1,
           legend.title = "Malignant state group",
           palette = "Set1", surv.scale = "percent", ggtheme = theme_bw()) #+
