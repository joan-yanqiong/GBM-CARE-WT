#################################
## Title: Spitzer Figure 4 panel f - comparison between longitudinal proportion difference of the MES state in CARE and Wang et al
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate a boxplot comparing the longitudinal proportion difference of the MES state in CARE and Wang et al
#################################

d <- mdata %>%
  filter(Patient %in% mgmt_md$Patient, Timepoint != "T3") %>%
  group_by(Sample, Patient, Timepoint, PvsR, State) %>%
  summarise(n = n()) %>%
  group_by(Sample, Patient, Timepoint, PvsR) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  ungroup() %>%
  left_join(mgmt_md %>%
              dplyr::select(Sample, MGMT, MGMT_exp_level),
            by = "Sample")

d_stats_summary <- acast(data = d, formula = Sample ~ State, value.var = "Freq")
d_stats_summary[is.na(d_stats_summary)] <- 0
d_stats_summary <- melt(d_stats_summary) %>%
  as_tibble()
colnames(d_stats_summary) <- c("Sample", "State", "Freq")

d_stats_summary <- d_stats_summary %>%
  left_join(d %>%
              ungroup() %>%
              dplyr::select(Sample, Patient, Timepoint, PvsR, MGMT, MGMT_exp_level) %>%
              filter(!duplicated(Sample)),
            by = "Sample")

d_stats_summary <- d_stats_summary %>%
  group_by(Patient, State) %>%
  mutate(MGMT_exp_level = MGMT_exp_level[Timepoint == "T1"]) %>%
  ungroup()

d_stats_summary_t2_vs_t1 <- d_stats_summary %>%
  group_by(State, Patient, MGMT_exp_level) %>%
  summarise(Diff = Freq[Timepoint == "T2"] - Freq[Timepoint == "T1"], .groups = "drop")

diaz_mes <- read.delim(paste0(DATA_ROOT, "diaz_longitudinal_mes_mgmt_t1_status.txt"), header = T, stringsAsFactors = F) %>%
  as_tibble()

diaz_mes$MGMT_T1 <- ifelse(diaz_mes$MGMT_T1 == "UM", "High", "Low")

diaz_mes <- diaz_mes %>%
  dplyr::select(UCSF_Pair, MGMT_T1, dMES) %>%
  dplyr::rename(Patient = UCSF_Pair, MGMT_exp_level = MGMT_T1, Diff = dMES) %>%
  mutate(Dataset = "Wang")

d <- d_stats_summary_t2_vs_t1 %>%
  filter(State == "MES") %>%
  dplyr::select(-State) %>%
  mutate(Dataset = "CARE")

d <- rbind(d, diaz_mes)

d <- d %>%
  filter(MGMT_exp_level != "Intermediate")

d_n <- d %>%
  group_by(MGMT_exp_level, Dataset) %>%
  summarise(n = n(), .groups = "drop")

d %>%
  ggboxplot(x = "MGMT_exp_level", y = "Diff", add = "jitter", color = "MGMT_exp_level") +
  facet_wrap(~Dataset, nrow = 1, scales = "free_x") +
  geom_text(data = d_n, mapping = aes(x = MGMT_exp_level, y = .3, label = paste0("n=", n)),
            size = 5, color = "black") +
  scale_color_brewer(name = "MGMT status", palette = "Set1") +
  scale_y_continuous(labels = percent, limits = c(-.2, .3), oob = squish) +
  xlab("Dataset") +
  ylab("Proportion") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "top", axis.text.x = element_text(size = 16)) +
  theme(panel.grid.major = element_line())

d %>%
  group_by(Dataset) %>%
  summarise(p = wilcox.test(Diff[MGMT_exp_level == "Low"], Diff[MGMT_exp_level == "High"])$p.value)
