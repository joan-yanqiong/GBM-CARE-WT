#################################
## Title: Spitzer Figure 2 panel d and extended data figure 3 panel e - conserved features across time points
## Date: 2025.02.24
## Author: Avishays Spitzer
## Description: Generate a heatmap and a graph showing the conserved features across time points
#################################

####################################################################################################################################
# Setup data
####################################################################################################################################

d <- feature_scores_all_per_sample %>%
  mutate(ID = paste0(Patient, Timepoint))

dm <- melt(d %>%
             select(ID, C1, C3, C5))
dm <- acast(dm, formula = ID ~ variable, value.var = "value")
colnames(dm) <- c("SCP_ECM", "SCP_Neuronal", "SCP_Glial")

m1 <- cbind(dt_malignant_composition_all_samples$T1_ABS,
            dt_overall_composition_all_samples$T1_ABS,
            dt_cell_cycle_all_samples$T1_ABS)
m2 <- cbind(dt_malignant_composition_all_samples$T2_ABS,
            dt_overall_composition_all_samples$T2_ABS,
            dt_cell_cycle_all_samples$T2_ABS)

m1 <- m1[, colnames(m1) %ni% c("T1_ABS_TME_Unresolved", "T1_ABS_TME_Other normal", "T1_ABS_TME_Other neuron", "T1_ABS_State_Unresolved")]
m2 <- m2[, colnames(m2) %ni% c("T2_ABS_TME_Unresolved", "T2_ABS_TME_Other normal", "T2_ABS_TME_Other neuron", "T2_ABS_State_Unresolved")]

rownames(m1) <- paste0(rownames(m1), "T1")
rownames(m2) <- paste0(rownames(m2), "T2")

m1 <- cbind(m1, dm[rownames(m1), ])
m2 <- cbind(m2, dm[rownames(m2), ])

d <- mdata %>%
  mutate(PvsR = gsub(" ", "", PvsR)) %>%
  mutate(PvsR = ifelse(PvsR == "Primary", "Primary", "Recurrent"))

valid_pts_p <- d %>%
  filter(PvsR == "Primary") %>%
  group_by(Patient) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(Patient)
valid_pts_p <- paste0(valid_pts_p, "T1")

valid_pts_r <- d %>%
  filter(PvsR == "Recurrent") %>%
  group_by(Patient) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(Patient)
valid_pts_r <- paste0(valid_pts_r, "T2")

m1 <- m1[valid_pts_p, ]
m2 <- m2[valid_pts_r, ]

cor_m1 <- cor(m1)
cor_m2 <- cor(m2)

hc1 <- fastcluster::hclust(d = as.dist(1 - cor_m1), method = "average")
hc2 <- fastcluster::hclust(d = as.dist(1 - cor_m2), method = "average")

t1_ft_cor <- lapply(colnames(m1), function(ft1) {
  res <- lapply(colnames(m1), function(ft2) {
    cr <- cor.test(m1[, ft1], m1[, ft2])
    tibble(Feature1 = ft1, Feature2 = ft2, cor = cr$estimate, pval = cr$p.value)
  })
  res <- do.call(rbind, res)
  return(res)
})
t1_ft_cor <- do.call(rbind, t1_ft_cor)

t1_ft_cor$Edge <- sapply(1:nrow(t1_ft_cor), function(i) {
  x <- sort(t1_ft_cor[i, c("Feature1", "Feature2")])
  paste0(x[1], "=", x[2])
})

t1_ft_cor <- t1_ft_cor %>%
  filter(Feature1 != Feature2)
t1_ft_cor <- t1_ft_cor %>%
  filter(!duplicated(Edge))

t1_ft_cor$padj <- p.adjust(t1_ft_cor$pval, "fdr")

t1_ft_cor$Edge <- gsub("T1_", "", t1_ft_cor$Edge)

t2_ft_cor <- lapply(colnames(m2), function(ft1) {
  res <- lapply(colnames(m2), function(ft2) {
    cr <- cor.test(m2[, ft1], m2[, ft2])
    tibble(Feature1 = ft1, Feature2 = ft2, cor = cr$estimate, pval = cr$p.value)
  })
  res <- do.call(rbind, res)
  return(res)
})
t2_ft_cor <- do.call(rbind, t2_ft_cor)

t2_ft_cor$Edge <- sapply(1:nrow(t2_ft_cor), function(i) {
  x <- sort(t2_ft_cor[i, c("Feature1", "Feature2")])
  paste0(x[1], "=", x[2])
})

t2_ft_cor <- t2_ft_cor %>%
  filter(Feature1 != Feature2)
t2_ft_cor <- t2_ft_cor %>%
  filter(!duplicated(Edge))

t2_ft_cor$padj <- p.adjust(t2_ft_cor$pval, "fdr")

t2_ft_cor$Edge <- gsub("T2_", "", t2_ft_cor$Edge)

t1_ft_cor <- t1_ft_cor %>%
  filter(pval < .05)

t2_ft_cor <- t2_ft_cor %>%
  filter(pval < .05)

####################################################################################################################################
# Figure EDF3e - Graph of longitudinally conserved features
####################################################################################################################################

int_fts <- intersect(t1_ft_cor$Edge[t1_ft_cor$pval < .05], t2_ft_cor$Edge[t2_ft_cor$pval < .05])

int_res <- t1_ft_cor %>%
  filter(Edge %in% int_fts) %>%
  left_join(t2_ft_cor %>%
              filter(Edge %in% int_fts),
            by = c("Edge"))

edges <- lapply(int_res$Edge[int_res$cor.x > 0], function(h) {
  
  h <- strsplit(x = as.character(h), split = "=")
  
  vs1 <- h[[1]][1]
  vs2 <- h[[1]][2]
  data.frame(S1 = vs1, S2 = vs2)
})
edges <- do.call(rbind, edges)

int_graph <- graph_from_data_frame(edges, directed = F)

colors <- RColorBrewer::brewer.pal(3, "Set1")

V(int_graph)$vertex.color[grep("SCP", rownames(int_graph[]))] <- colors[1]
V(int_graph)$vertex.color[grep("TME", rownames(int_graph[]))] <- colors[2]
V(int_graph)$vertex.color[grep("State", rownames(int_graph[]))] <- colors[3]

V(int_graph)$label <- sapply(strsplit(rownames(int_graph[]), split = "_"), function(x) x[[length(x)]][1])

plot(int_graph,
     vertex.color = V(int_graph)$vertex.color, vertex.label.color="black",
     vertex.label = V(int_graph)$label,
     edge.label.cex = 1,
     edge.curved = F,
     vertex.size = 10, vertex.label.cex = 1.5)

####################################################################################################################################
# Figure 2d - Combined feature correlation matrix
####################################################################################################################################

comb_cor_matrix <- psych::lowerUpper(cor_m1, cor_m2)

dm <- melt(comb_cor_matrix) %>%
  as_tibble() %>%
  mutate(Var1 = factor(as.character(Var1), hc2$labels[hc2$order]),
         Var2 = factor(as.character(Var2), hc1$labels[hc2$order]))

dm$E1 <- gsub("T2_", "", as.character(dm$Var1))
dm$E2 <- gsub("T1_", "", as.character(dm$Var2))

fts <- unique(c(edges$S1, edges$S2))

dm <- dm %>%
  filter(E1 %in% fts & E2 %in% fts)

ggplot(dm, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(name = "Correlation", low = "dodgerblue", mid = "white", high = "red",
                       limits = c(-.5, .5), oob = squish, labels = c("-.5", "", "0", "", ".5")) +
  xlab("T2 features") +
  ylab("T1 features") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank(),
                 axis.text.x = element_text(size = 16, angle = 90), axis.text.y = element_text(size = 16))
