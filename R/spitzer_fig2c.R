#################################
## Title: Spitzer Figure 2 panel c - survival curves
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate survival curves for the LMF-GN composition group and MES/Hypoxia malignant state group
#################################

####################################################################################################################################
# LMF - GN
####################################################################################################################################

d <- gbm_subtypes_tbl %>%
  
  d$isGN <- d$CompCluster == "LP - GN"
  
  km_fit <- survfit(Surv(SurgicalInterval2, VT) ~ isGN, data = d)
  ggsurvplot(fit = km_fit, data = d,
             pval = T,
             risk.table = T,
             xlim = c(0, 30),
             break.x.by = 10,    
             xlab = "Time in months",
             size = 1,
             legend.title = "GN vs. not GN",
             palette = "Set1", surv.scale = "percent", ggtheme = theme_bw()) #+
  
  ####################################################################################################################################
  # MES/Hypoxia
  ####################################################################################################################################
  
  d <- gbm_subtypes_tbl %>%
    filter(Timepoint == "T2") %>%
    mutate(MalCluster = as.character(MalCluster))
  
  d$MalCluster[is.na(d$MalCluster)] <- "Unknown"
  
  d$isMES <- d$MalCluster == "MES/Hypoxia"
  
  km_fit <- survfit(Surv(SurgicalInterval2, VT) ~ isMES, data = d)
  ggsurvplot(fit = km_fit, data = d,
             pval = T,
             risk.table = T,
             xlim = c(0, 30),
             break.x.by = 10,    
             xlab = "Time in months",
             size = 1,
             legend.title = "MES/Hypoxia vs. not MES/Hypoxia",
             palette = "Set1", surv.scale = "percent", ggtheme = theme_bw()) #+
  