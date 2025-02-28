#################################
## Title: Extended Data Figure 1 panel c in Nomura et al - quality control metrics
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Sample-level copy-number matrix (per chromosome)
#################################

d <- cna_sig_per_chr %>%
  filter(CellID %in% mdata$CellID)

d <- d %>%
  group_by(Patient, Timepoint, Sample, CHR) %>%
  summarise(Sig = mean(Sig))
d$ID <- paste0(d$Patient, d$Timepoint)

d$Sig[abs(d$Sig) < .05] <- 0

m <- acast(d, ID ~ CHR, value.var = "Sig")

hc <- fastcluster::hclust(d = as.dist(1 - cor(t(m))), method = "average")

d$ID <- factor(as.character(d$ID), hc$labels[hc$order])

ggplot(d, aes(x = CHR, y = ID, fill = Sig)) +
  facet_grid(cols = vars(CHR), scales = "free_x") +
  geom_tile() +
  scale_fill_gradient2(name = "", low = "dodgerblue", mid = "white", high = "red",
                       labels = c("DEL", "", "", "", "AMP"), limits = c(-.2, .2), oob = squish) +
  guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 1, ticks.colour = "black")) +
  xlab("Chromosomes") +
  ylab("Samples") +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90),
                 axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                 legend.text = element_text(size = 16),
                 strip.text = element_text(size = 10))
