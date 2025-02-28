#################################
## Title: Figure 2 panel f in Nomura et al - example of hybrid states (AC, GPC)
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a gene-expression heatmap of hybrid cells
#################################

####################################################################################################################################
# Select the singular cells
####################################################################################################################################

singular_states <- state_data %>%
  ungroup() %>%
  filter(Program %in% c("MP_4_AC", "MP_8_GPC"), p.sig == T) %>%
  filter(CellID %ni% state_data_hybrid$CellID)

singular_states <- singular_states %>%
  group_by(Program) %>%
  arrange(Program, desc(Score)) %>%
  top_n(100, Score) %>%
  ungroup()

####################################################################################################################################
# Select the AC-GPC hybrids
####################################################################################################################################

ac_gpc_hybrids <- state_data %>%
  ungroup() %>%
  filter(Program %in% c("MP_4_AC", "MP_8_GPC"), p.sig == T) %>%
  group_by(CellID) %>%
  filter("MP_4_AC" %in% Program & "MP_8_GPC" %in% Program) %>%
  ungroup() %>%
  filter(CellID %in% state_data_hybrid$CellID[state_data_hybrid$N == 2])

ac_gpc_hybrids <- ac_gpc_hybrids %>%
  filter(Program == "MP_4_AC") %>%
  rename(AC = Score) %>%
  select(CellID, AC) %>%
  left_join(ac_gpc_hybrids %>%
              filter(Program == "MP_8_GPC") %>%
              rename(GPC = Score) %>%
              select(CellID, GPC),
            by = "CellID")
ac_gpc_hybrids$Diff <- abs(ac_gpc_hybrids$AC - ac_gpc_hybrids$GPC)

ac_gpc_hybrids <- ac_gpc_hybrids %>%
  mutate(Sum = AC + GPC)

ac_gpc_hybrids <- ac_gpc_hybrids %>%
  arrange(Diff)

gghistogram(data = ac_gpc_hybrids, x = "Diff", bins = 50)

ac_gpc_hybrids <- ac_gpc_hybrids %>%
  filter(Diff < .01) %>%
  arrange(desc(Sum)) %>%
  top_n(100, Diff)

####################################################################################################################################
# Generate the heatmap
####################################################################################################################################

hybrid_vs_singular_cells <- rbind(singular_states %>%
                                    select(CellID, Program) %>%
                                    mutate(State = case_when(Program == "MP_4_AC" ~ "AC",
                                                             Program == "MP_8_GPC" ~ "GPC")) %>%
                                    select(CellID, State),
                                  ac_gpc_hybrids %>%
                                    select(CellID) %>%
                                    mutate(State = "AC-GPC"))

hybrid_vs_singular_cells$Sample <- get_sname(hybrid_vs_singular_cells$CellID, split = "-")

mp_genes <- rbind(tibble(Gene = MP_list_named$MP_8_GPC, MP = "GPC"),
                  tibble(Gene = MP_list_named$MP_4_AC, MP = "AC"))

m <- lapply(unique(hybrid_vs_singular_cells$Sample), function(sname) {
  print(sname)
  res <- umi_data_all[[sname]][mp_genes$Gene, hybrid_vs_singular_cells$CellID[hybrid_vs_singular_cells$Sample == sname]]
  if(is.null(dim(res))) {
    res <- as.matrix(res)
    colnames(res) <- hybrid_vs_singular_cells$CellID[hybrid_vs_singular_cells$Sample == sname]
    rownames(res) <- mp_genes$Gene
  }
  return(res)
})
m <- do.call(cbind, m)

m <- umi2upm(m)
m <- log2(m / 10 + 1)
m <- rowcenter(m)

dm <- melt(m) %>%
  as_tibble()
colnames(dm) <- c("Gene", "CellID", "value")

dm <- dm %>%
  left_join(hybrid_vs_singular_cells %>%
              select(CellID, State),
            by = "CellID") %>%
  mutate(State = factor(State, c("AC", "AC-GPC", "GPC")))

dm <- dm %>%
  left_join(mp_genes %>%
              select(Gene, MP),
            by = "Gene") %>%
  mutate(MP = factor(MP, c("AC", "GPC")))

ord1 <- ac_gpc_hybrids %>%
  mutate(Diff = GPC - AC) %>%
  arrange(Diff) %>%
  pull(CellID)

ord2 <- singular_states %>%
  pull(CellID)

dm$CellID <- factor(dm$CellID, c(ord1, ord2))

ggplot(dm, aes(x = CellID, y = Gene, fill = value)) +
  facet_grid(rows = vars(MP), cols = vars(State), scales = "free") +
  geom_tile() +
  scale_fill_distiller(name = "Centered log2-expression", palette = "RdBu", direction = -1,
                       breaks = seq(-3, 3),
                       limits = c(-3, 3), oob = squish, guide = "colorsteps") +
  xlab("Cells") +
  ylab("Genes") +
  guides(fill = guide_colorsteps(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 7))
