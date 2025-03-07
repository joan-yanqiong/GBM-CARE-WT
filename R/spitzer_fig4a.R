#################################
## Title: Spitzer Figure 4 panel a - MGMT expression across time points and MGMT status
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate a boxplot of MGMT expression across time points and MGMT status
#################################

####################################################################################################################################
# Setup
####################################################################################################################################

valid_samples <- sample_data %>%
  filter(IDHstatus == "WT", Timepoint != "T3", Patient %in% pt_pairs$Patient) %>%
  group_by(Patient) %>%
  filter("Primary" %in% PvsR & "1st recurrence" %in% PvsR)

mgmt_status_tbl <- read.csv(file = paste0(TABLES_ROOT, "Metadata sheet P-R GBM 03272023 - Clinical info.csv"), header = T, stringsAsFactors = F) %>%
  as_tibble() %>%
  select(Sample.ID.Avishay.named, MGMT) %>%
  rename(ID = Sample.ID.Avishay.named) %>%
  filter(ID %in% valid_samples$ID) %>%
  mutate(MGMT = case_when(is.na(MGMT) ~ "NOS",
                          MGMT %in% c("Na", "UA") ~ "NOS",
                          MGMT == "Methylated" ~ "MET",
                          MGMT == "Unmethylated" ~ "UM"))

mgmt_status_vec <- setNames(mgmt_status_tbl$MGMT, mgmt_status_tbl$ID)

mgmt_md <- gbm_subtypes_tbl %>%
  filter(Sample %in% valid_samples$Sample)

mgmt_md$MGMT <- mgmt_status_vec[as.character(mgmt_md$ID)]

mgmt_md$MGMT_orig <- mgmt_md$MGMT

mgmt_md <- mgmt_md %>%
  group_by(Patient) %>%
  mutate(MGMT = MGMT[Timepoint == "T1"]) %>%
  ungroup()

# Compute MGMT expression level per sample
mgmt_exp_vec <- sapply(mgmt_md$Sample, function(sname) {
  print(sname)
  cellids <- mdata %>%
    filter(Sample == sname) %>%
    pull(CellID)
  exp <- umi_data_all[[sname]][, cellids]
  exp <- umi2upm(exp)
  log2(mean(exp["MGMT", ]) + 1)
})

mgmt_md$MGMT_exp <- mgmt_exp_vec[mgmt_md$Sample]

####################################################################################################################################
# Figure 4a - MGMT expression across time points and MGMT status
####################################################################################################################################

mgmt_md %>%
  filter(MGMT != "NOS") %>%
  mutate(MGMT = factor(MGMT, c("MET", "UM", "NOS"))) %>%
  ggplot(aes(x = Timepoint, y = MGMT_exp)) +
  facet_wrap(~MGMT, scales = "free_x") +
  geom_boxplot(aes(color = Timepoint), size = 1) +
  geom_line(aes(group = Patient), linetype = "dashed", color = "black", size = .5) +
  xlab("Timepoint") +
  ylab("Mean MGMT expression\n(per sample)") +
  scale_color_manual(values = timepoint_color_vec) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line(), legend.position = "none")

####################################################################################################################################
# Stats
####################################################################################################################################

d <- mgmt_md %>%
  filter(MGMT != "NOS")

wilcox.test(x = mgmt_md %>% filter(MGMT == "MET", Timepoint == "T1") %>% pull(MGMT_exp),
            y = mgmt_md %>% filter(MGMT == "MET", Timepoint == "T2") %>% pull(MGMT_exp))

wilcox.test(x = mgmt_md %>% filter(MGMT == "UM", Timepoint == "T1") %>% pull(MGMT_exp),
            y = mgmt_md %>% filter(MGMT == "UM", Timepoint == "T2") %>% pull(MGMT_exp))

wilcox.test(x = mgmt_md %>% filter(MGMT == "MET", Timepoint == "T1") %>% pull(MGMT_exp),
            y = mgmt_md %>% filter(MGMT == "UM", Timepoint == "T1") %>% pull(MGMT_exp))

wilcox.test(x = mgmt_md %>% filter(MGMT == "MET", Timepoint == "T2") %>% pull(MGMT_exp),
            y = mgmt_md %>% filter(MGMT == "UM", Timepoint == "T2") %>% pull(MGMT_exp))

wilcox.test(x = mgmt_md %>% filter(MGMT == "MET") %>% pull(MGMT_exp),
            y = mgmt_md %>% filter(MGMT == "UM") %>% pull(MGMT_exp))
