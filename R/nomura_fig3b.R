#################################
## Title: Figure 3 panel b in Nomura et al - Inter vs. intra tumor heterogeneity
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: compare inter vs. intra-tumor variability (i.e. variability of meta-programs and baseline profiles across tumor samples)
#################################

mp_list <- c(MP_list_named,
             list(BP_ECM = pca_gene_sigs$C1,
                  BP_Neuronal = pca_gene_sigs$C3,
                  BP_Glial = pca_gene_sigs$C5))

mp_vs_bp_scores <- lapply(names(mp_list), function(mpname) {
  print(mpname)
  mp <- mp_list[[mpname]]
  res <- lapply(state_pb_profiles, function(pb) {
    pb <- pb[rownames(pb) %in% mp, ]
    tibble(ID = colnames(pb),
           Sample = sapply(strsplit(ID, split = "_"), function(x) x[[1]][1]),
           State = sapply(strsplit(ID, split = "_"), function(x) x[[2]]),
           MP = mpname, Score = colMeans(pb))
  })
  res <- do.call(rbind, res)
  return(res)
})
mp_vs_bp_scores <- do.call(rbind, mp_vs_bp_scores)

x <- mp_vs_bp_scores %>%
  group_by(Sample, MP) %>%
  summarise(SD = sd(Score), .groups = "drop") %>%
  filter(!is.na(SD))

x <- x %>%
  group_by(MP) %>%
  summarise(Intra = mean(SD), .groups = "drop")

y <- mp_vs_bp_scores %>%
  group_by(State, MP) %>%
  summarise(SD = sd(Score), .groups = "drop")

y <- y %>%
  group_by(MP) %>%
  summarise(Inter = mean(SD), .groups = "drop")

d <- x %>%
  left_join(y, by = "MP")

d <- d %>%
  filter(MP %in% c("MP_2_OPC", "MP_4_AC", "MP_5_Hypoxia", "MP_6_MES", "MP_7_NPC", "MP_8_GPC", "MP_9_ExN",
                   "BP_ECM", "BP_Neuronal", "BP_Glial"))

d$Type <- sapply(strsplit(d$MP, split = "_"), function(x) x[[1]][1])

ggplot(d, aes(x = Intra, y = Inter)) +
  geom_point(aes(fill = Type), color = "black", shape = 21, size = 8) +
  geom_text(aes(label = MP), size = 8) +
  scale_fill_brewer(name = "", palette = "Set1", labels = c("BP" = "Baseline profile", "MP" = "Meta-program")) +
  scale_x_continuous(limits = c(.5, 1)) +
  scale_y_continuous(limits = c(.5, 1.3)) +
  xlab("Intra-tumor variability") +
  ylab("Inter-tumor variability") +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

t.test(d$Inter[d$Type == "BP"], d$Inter[d$Type == "MP"], var.equal = T, alternative = "greater")
t.test(d$Intra[d$Type == "BP"], d$Intra[d$Type == "MP"], var.equal = T, alternative = "less")
