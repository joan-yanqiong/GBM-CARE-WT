#################################
## Title: Extended Data Figure 3 panel b-d in Nomura et al - comparison with Neftel data
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Compare CARE states with Neftel states
#################################

####################################################################################################################################
# Prepare data for plotting
####################################################################################################################################

nl_exp_data <- read.delim("IDHwtGBM.processed.SS2.logTPM.txt", header = T)

rownames(nl_exp_data) <- nl_exp_data[, 1]
nl_exp_data <- nl_exp_data[, -1]
nl_exp_data <- as.matrix(nl_exp_data)
colnames(nl_exp_data) <- gsub("\\.", "-", colnames(nl_exp_data))

nl_meta_data <- read.delim("IDHwt.GBM.Metadata.SS2.txt", header = T, stringsAsFactors = F)

nl_meta_data <- nl_meta_data[-1, ]
nl_meta_data <- as_tibble(nl_meta_data)

nl_meta_data <- nl_meta_data %>%
  filter(CellAssignment == "Malignant")

sigs <- c(setNames(Signatures_GBM, paste0("NL_", names(Signatures_GBM))),
          MP_list_named)

nl_scores <- lapply(unique(nl_meta_data$Sample), function(s_name) {
  
  d <- nl_meta_data %>% filter(Sample == s_name)
  
  x <- nl_exp_data[, d$NAME]
  dim(x)
  
  res <- scalop::sigScores(m = x, sigs = sigs, conserved.genes = .25)
  res <- as_tibble(res, rownames = "NAME")
  
  res$Sample <- s_name
  
  return(res)
})
nl_scores <- do.call(rbind, nl_scores)

nl_states <- c("NL_AC", "NL_MES1", "NL_MES2", "NL_OPC", "NL_NPC1", "NL_NPC2")
mp_states <- c("MP_2_OPC", "MP_4_AC", "MP_5_Hypoxia", "MP_6_MES", "MP_7_NPC", "MP_8_GPC", "MP_9_ExN", "MP_13_Cilia")

nl_scores$MaxNLscore <- apply(nl_scores[, nl_states], 1, max)
nl_scores$MaxNLstate <- apply(nl_scores[, nl_states], 1, function(x) nl_states[which.max(x)])

nl_scores$MaxMPscore <- apply(nl_scores[, mp_states], 1, max)
nl_scores$MaxMPstate <- apply(nl_scores[, mp_states], 1, function(x) mp_states[which.max(x)])

nl_scores_care <- score_within_samples(umi_data_list = umi_data_all, md = mdata %>% filter(Sample %in% pt_pairs$Sample), sigs = setNames(Signatures_GBM, paste0("NL_", names(Signatures_GBM))))

nl_scores_care <- nl_scores_care %>%
  left_join(MP_scores %>%
              select(CellID, starts_with("MP_")), by = "CellID")

nl_scores_care$MaxNLscore <- apply(nl_scores_care[, nl_states], 1, max)
nl_scores_care$MaxNLstate <- apply(nl_scores_care[, nl_states], 1, function(x) nl_states[which.max(x)])

nl_scores_care$MaxMPscore <- apply(nl_scores_care[, mp_states], 1, max)
nl_scores_care$MaxMPstate <- apply(nl_scores_care[, mp_states], 1, function(x) mp_states[which.max(x)])

####################################################################################################################################
# Figure EDF3b - NL vs. CARE scores correlation
####################################################################################################################################

cor_m <- cor(nl_scores_care %>%
               select(starts_with("NL_")),
             nl_scores_care %>%
               select(all_of(mps)))

hc1 <- fastcluster::hclust(d = dist(cor_m, "euclidean"))
hc2 <- fastcluster::hclust(d = dist(t(cor_m), "euclidean"))

cor_m <- melt(cor_m) %>%
  as_tibble()

cor_m$Var1 <- factor(as.character(cor_m$Var1), hc1$labels[hc1$order])
cor_m$Var2 <- factor(as.character(cor_m$Var2), hc2$labels[hc2$order])

ggplot(cor_m, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(name = "Correlation", low = "dodgerblue", mid = "white", high = "red", limits = c(-.75, .75), oob = squish,
                       breaks = c(-.75, -.5, -.25, 0, .25, .5, .75), labels = c("-.75", "", "", "0", "", "", ".75")) +
  xlab("CARE MPs") +
  ylab("Neftel MPs") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 20, angle = 90))

####################################################################################################################################
# Figure EDF3c - NL vs. CARE Jaccard
####################################################################################################################################

mps <- c("MP_2_OPC", "MP_3_CC", "MP_4_AC", "MP_5_Hypoxia", "MP_6_MES", "MP_7_NPC", "MP_8_GPC", "MP_9_ExN", "MP_10_Stress1", "MP_13_Cilia", "MP_14_NRGN", "MP_15_Stress2")

jac_m <- jaccard(MP_list_named[mps], Signatures_GBM)

hc1 <- fastcluster::hclust(d = dist(jac_m, "euclidean"))
hc2 <- fastcluster::hclust(d = dist(t(jac_m), "euclidean"))

jac_m <- melt(jac_m) %>%
  as_tibble()

jac_m$Var1 <- factor(as.character(jac_m$Var1), hc1$labels[hc1$order])
jac_m$Var2 <- factor(as.character(jac_m$Var2), hc2$labels[hc2$order])

ggplot(jac_m, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(name = "Jaccard", low = "white", high = "red", limits = c(0, .15), oob = squish) +
  xlab("CARE MPs") +
  ylab("Neftel MPs") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 20, angle = 90))

####################################################################################################################################
# Figure EDF3d - Compare NL scores with CARE scores
####################################################################################################################################

res1 <- rbind(nl_scores %>%
                filter(MaxNLstate == "NL_AC" & MaxMPstate == "MP_4_AC") %>%
                mutate(State = "AC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores %>%
                filter(MaxNLstate == "NL_MES1" & MaxMPstate == "MP_6_MES") %>%
                mutate(State = "MES") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores %>%
                filter(MaxNLstate == "NL_MES2" & MaxMPstate == "MP_5_Hypoxia") %>%
                mutate(State = "Hypoxia") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores %>%
                filter(MaxNLstate == "NL_OPC" & MaxMPstate == "MP_2_OPC") %>%
                mutate(State = "OPC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores %>%
                filter(MaxNLstate == "NL_NPC2" & MaxMPstate == "MP_7_NPC") %>%
                mutate(State = "NPC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt()) %>%
  as_tibble() %>%
  mutate(Dataset = "Neftel")

res2 <- rbind(nl_scores_care %>%
                filter(MaxNLstate == "NL_AC" & MaxMPstate == "MP_4_AC") %>%
                mutate(State = "AC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores_care %>%
                filter(MaxNLstate == "NL_MES1" & MaxMPstate == "MP_6_MES") %>%
                mutate(State = "MES") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores_care %>%
                filter(MaxNLstate == "NL_MES2" & MaxMPstate == "MP_5_Hypoxia") %>%
                mutate(State = "Hypoxia") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores_care %>%
                filter(MaxNLstate == "NL_OPC" & MaxMPstate == "MP_2_OPC") %>%
                mutate(State = "OPC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt(),
              nl_scores_care %>%
                filter(MaxNLstate == "NL_NPC2" & MaxMPstate == "MP_7_NPC") %>%
                mutate(State = "NPC") %>%
                select(State, MaxNLscore, MaxMPscore) %>%
                rename(Neftel = MaxNLscore, CARE = MaxMPscore) %>%
                melt()) %>%
  as_tibble() %>%
  mutate(Dataset = "CARE")

res <- rbind(res1, res2)

res_stats <- res %>%
  group_by(Dataset, State, variable) %>%
  summarise(Mean = mean(value), .groups = "drop")

p1 <- res %>%
  filter(Dataset == "CARE") %>%
  ggplot(aes(x = value, y = after_stat(ndensity), fill = variable)) +
  geom_density(color = "black", alpha = .3) +
  geom_vline(data = res_stats %>%
               filter(Dataset == "CARE"),
             mapping = aes(xintercept = Mean, color = variable), linetype = "dashed", size = 1, show.legend = F) +
  facet_grid(rows = vars(Dataset), cols = vars(State)) +
  scale_fill_discrete(name = "Signature origin") +
  xlab("Score") +
  ylab("Density (scaled to 1)") +
  theme_gbm_pvsr()

p2 <- res %>%
  filter(Dataset == "Neftel") %>%
  ggplot(aes(x = value, y = after_stat(ndensity), fill = variable)) +
  geom_density(color = "black", alpha = .3) +
  geom_vline(data = res_stats %>%
               filter(Dataset == "Neftel"),
             mapping = aes(xintercept = Mean, color = variable), linetype = "dashed", size = 1, show.legend = F) +
  facet_grid(rows = vars(Dataset), cols = vars(State)) +
  scale_fill_discrete(name = "Signature origin") +
  xlab("Score") +
  ylab("Density (scaled to 1)") +
  theme_gbm_pvsr()

p1 + p2 + plot_layout(nrow = 2)
