#################################
## Title: Extended Data Figure 4 panel a in Nomura et al - fraction of cycling cells across states
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a boxplot showing the fraction of cycling cells across states
#################################

d <- mdata %>%
  filter(State %in% c("Cilia", "AC", "MES", "Hypoxia", "Stress", "GPC", "OPC", "NPC", "Neuron"))
table(d$State)
table(d$State, d$isCC)
table(d$State, d$isCC) / rowSums(table(d$State, d$isCC))

d$PvsR <- gsub(" ", "", d$PvsR)
d$PvsR[grep("Recurrent", d$PvsR)] <- "Recurrence"
table(d$PvsR)

# Exclude samples with less than 50 malignant cells (to reduce noise)
d <- d %>%
  group_by(Sample, PvsR, State) %>%
  summarise(n = sum(isCC), N = n(), Freq = n / N, .groups = "drop") %>%
  filter(N >= 50)

d %>%
  group_by(State) %>%
  summarise(Mean = mean(Freq), SD = sd(Freq), n = n(), SE = SD / sqrt(n))

d <- d %>%
  mutate(State_fct = case_when(State == "Cilia" ~ "Cilia-like",
                               State == "AC" ~ "AC-like",
                               State == "MES" ~ "MES-like",
                               State == "Hypoxia" ~ "Hypoxia",
                               State == "Stress" ~ "Stress",
                               State == "GPC" ~ "GPC-like",
                               State == "OPC" ~ "OPC-like",
                               State == "NPC" ~ "NPC-like",
                               State == "Neuron" ~ "Neuron-like")) %>%
  mutate(State_fct = factor(State_fct, c("Cilia-like", "AC-like", "MES-like", "Hypoxia",
                                         "Stress", "GPC-like", "OPC-like", "NPC-like",
                                         "Neuron-like")))

d_ord <- d %>%
  group_by(State_fct) %>%
  summarise(Median = median(Freq), Mean = mean(Freq), .groups = "drop") %>%
  arrange(Median)

d$State_fct <- factor(as.character(d$State_fct), as.character(d_ord$State_fct))

ggboxplot(data = d, x = "State_fct", y = "Freq", add = "jitter", xlab = "State", ylab = "% cycling cells") +
  scale_y_continuous(labels = percent, limits = c(0, .4), oob = squish) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())
