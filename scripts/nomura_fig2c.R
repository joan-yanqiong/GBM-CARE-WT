#################################
## Title: Figure 2 panel c in Nomura et al - gene-expression heatmap of new states
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a gene-expression heatmap of new states
#################################

new_states <- mdata %>%
  filter(State_full %in% c("MP_8_GPC", "MP_9_ExN", "MP_13_Cilia"))

new_states_scores <- state_data %>%
  ungroup() %>%
  filter(CellID %in% new_states$CellID, Program %in% c("MP_8_GPC", "MP_9_ExN", "MP_13_Cilia"), p.sig == T)

is_dup <- new_states_scores %>%
  group_by(CellID) %>%
  filter(duplicated(CellID)) %>%
  pull(CellID)

new_states_scores <- new_states_scores %>%
  filter(CellID %ni% is_dup)

new_states_scores$State <- new_states_scores$Class

selected_cells <- new_states_scores %>%
  group_by(Program) %>%
  arrange(desc(Score)) %>%
  top_frac(.05, Score) %>%
  ungroup()

reference_cells <- mdata %>%
  filter(State %in% c("AC", "MES", "Hypoxia", "OPC", "NPC")) %>%
  group_by(State) %>%
  sample_n(1000)

cells <- rbind(selected_cells %>%
                 select(CellID, Sample),
               reference_cells %>%
                 ungroup() %>%
                 select(CellID, Sample))

mp_genes <- rbind(tibble(Gene = MP_list_named$MP_13_Cilia, MP = "Cilia"),
                  tibble(Gene = MP_list_named$MP_8_GPC, MP = "GPC"),
                  tibble(Gene = MP_list_named$MP_9_ExN, MP = "Neuron"))

m <- lapply(unique(cells$Sample), function(sname) {
  print(sname)
  res <- umi_data_all[[sname]][mp_genes$Gene[mp_genes$Gene %in% rownames(umi_data_all[[sname]])], cells$CellID[cells$Sample == sname]]
  if(is.null(dim(res))) {
    res <- as.matrix(res)
    colnames(res) <- cells$CellID[cells$Sample == sname]
    rownames(res) <- mp_genes$Gene[mp_genes$Gene %in% rownames(umi_data_all[[sname]])]
  }
  return(res)
})
m <- do.call(cbind, m)

m <- umi2upm(m)
m <- log2(m / 10 + 1)
m <- rowcenter(m)

m <- m[, selected_cells$CellID]

dm <- melt(m) %>%
  as_tibble()
colnames(dm) <- c("Gene", "CellID", "value")

dm <- dm %>%
  left_join(selected_cells %>%
              select(CellID, State),
            by = "CellID") %>%
  mutate(State = factor(State, c("Cilia", "GPC", "Neuron")))

dm <- dm %>%
  left_join(mp_genes %>%
              select(Gene, MP),
            by = "Gene") %>%
  mutate(MP = factor(MP, c("Cilia", "GPC", "Neuron")))

ggplot(dm, aes(x = CellID, y = Gene, fill = value)) +
  facet_grid(rows = vars(MP), cols = vars(State), scales = "free") +
  geom_tile() +
  scale_fill_distiller(name = "Centered log2-expression", palette = "RdBu", direction = -1,
                       breaks = seq(-5, 5),
                       limits = c(-5, 5), oob = squish, guide = "colorsteps") +
  xlab("Cells") +
  ylab("Genes") +
  guides(fill = guide_colorsteps(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 7))
