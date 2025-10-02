#################################
## Title: Extended Data Figure 2h in Spitzer et al - Cox regression analysis
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce forrest plots of the cox regression analysis
#################################

####################################################################################################################################
# Composition group
####################################################################################################################################

d <- gbm_subtypes_tbl %>%
  filter(Timepoint == "T2")

d$Nrecurrences <- as.factor(d$Nrecurrences)
table(d$Nrecurrences)

cox_fit <- coxph(formula = Surv(SurgicalInterval2, VT) ~ CompCluster + Age + Gender + MGMT + Nrecurrences, data = d)

survminer::ggforest(model = cox_fit, data = d, fontsize = 2)

####################################################################################################################################
# Malignant state group
####################################################################################################################################

d <- gbm_subtypes_tbl %>%
  filter(Timepoint == "T2") %>%
  mutate(MalCluster = as.character(MalCluster))

d$MalCluster[is.na(d$MalCluster)] <- "Unknown"

cox_fit <- coxph(formula = Surv(SurgicalInterval2, VT) ~ MalCluster + Age + Gender + MGMT + Nrecurrences, data = d)
summary(cox_fit)

survminer::ggforest(model = cox_fit, data = d, fontsize = 2)
