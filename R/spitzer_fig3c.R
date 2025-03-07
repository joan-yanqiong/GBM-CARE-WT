#################################
## Title: Spitzer Figure 3 panel c - bbserved vs. expected conservation
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Compute the observed vs. expected conservation of transcriptional groups across timepoints
#################################

d_pres <- rbind(gbm_subtypes_tbl %>%
                  group_by(Patient) %>%
                  summarise(CompCluster = CompCluster[Timepoint == "T1"] == CompCluster[Timepoint == "T2"],
                            .groups = "drop") %>%
                  melt(id.vars = "Patient"),
                gbm_subtypes_tbl %>%
                  filter(!is.na(MalCluster)) %>%
                  group_by(Patient) %>%
                  filter("T1" %in% Timepoint & "T2" %in% Timepoint) %>%
                  summarise(MalCluster = MalCluster[Timepoint == "T1"] == MalCluster[Timepoint == "T2"],
                            .groups = "drop") %>%
                  melt(id.vars = "Patient"),
                gbm_subtypes_tbl %>%
                  filter(!is.na(SCP)) %>%
                  group_by(Patient) %>%
                  filter("T1" %in% Timepoint & "T2" %in% Timepoint) %>%
                  summarise(SCP = SCP[Timepoint == "T1"] == SCP[Timepoint == "T2"],
                            .groups = "drop") %>%
                  melt(id.vars = "Patient")) %>%
  as_tibble()

d_prop <- gbm_subtypes_tbl %>%
  summarise(CompCluster = 1 / length(levels(CompCluster)),
            MalCluster = 1 / length(levels(MalCluster)),
            SCP = 1 / length(levels(SCP)),
            # SCP2 = 1 / length(levels(SCP2)),
            .groups = "drop") %>%
  melt() %>%
  as_tibble()

d_pres_stats <- d_pres %>%
  group_by(variable) %>%
  summarise(x = sum(value, na.rm = T), n = n(), Freq = x / n, .groups = "drop")

d_pres_stats$Exp <- round(d_pres_stats$n * d_prop$value, 0)

binom.res <- lapply(1:nrow(d_pres_stats), function(i) binom.test(x = d_pres_stats$x[i], n = d_pres_stats$n[i], p = d_prop$value[i]))
binom.res <- lapply(binom.res, function(x) tibble(estimate = x$estimate, ci1 = x$conf.int[1], ci2 = x$conf.int[2], p.val = x$p.value))
binom.res <- do.call(rbind, binom.res)
binom.res$variable <- d_pres_stats$variable

binom.res$p.adj <- p.adjust(binom.res$p.val, "holm")

binom.res <- binom.res %>%
  mutate(plabel = case_when(p.adj < .001 ~ "***",
                            p.adj < .01 ~ "**",
                            p.adj < .05 ~ "*",
                            TRUE ~ ""))

ggplot(binom.res, aes(x = variable, y = estimate, fill = variable)) +
  geom_bar(color = "black", stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ci1, ymax = ci2), width = .2) +
  geom_text(aes(label = plabel, y = ci2 + .05), size = 8) +
  geom_point(data = d_prop, mapping = aes(x = variable, y = value), color = "red", fill = "black", shape = 21, size = 5, inherit.aes = F) +
  scale_fill_discrete(name = "") +
  scale_y_continuous(labels = percent) +
  scale_x_discrete(labels = c("SCP" = "Assigned\nSCP", "MalCluster" = "Malignant\nstate cluster", "CompCluster" = "Composition\ncluster")) +
  xlab("") +
  ylab("Observed conservation [%]") +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "none") +
  theme(panel.grid.major = element_line(), axis.text.x = element_text(size = 16))
