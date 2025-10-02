#################################
## Title: Extended Data Figure 5c in Spitzer et al - logistc regression for predicting MGMT methylation status from MGMT gene expression level
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a logistc regression model predicting MGMT methylation status from MGMT gene expression level
#################################

d <- mgmt_md %>%
  filter(MGMT != "NOS", Timepoint == "T1") %>%
  mutate(MGMT = case_when(MGMT == "MET" ~ 1,
                          TRUE ~ 0))

glm_fit <- glm(formula = MGMT ~ MGMT_exp, data = d, family = "binomial")
summary(glm_fit)

d <- mgmt_md

mgmt_md$MGMT_pred <- predict(glm_fit, newdata = d, type = "response")

ggscatter(data = mgmt_md, x = "MGMT_exp", y = "MGMT_pred", xlab = "Mean MGMT expression (per sample)", ylab = "Prabability of methylation", color = "MGMT") +
  geom_hline(yintercept = .75, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = .25, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = 3, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = mgmt_color_vec) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

mgmt_md <- mgmt_md %>%
  mutate(MGMT_exp_level = case_when(MGMT_pred > .75 ~ "Low",
                                    MGMT_pred < .25 ~ "High",
                                    TRUE ~ "Intermediate")) %>%
  mutate(MGMT_exp_level = as.factor(MGMT_exp_level))
