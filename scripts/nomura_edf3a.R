#################################
## Title: Extended Data Figure 3 panel a in Nomura et al - heatmap of meta-programs including low-quality/technical MPs
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a heatmap of meta-programs including the ones reflecting technical variability
#################################

limits <- lengths(malignant_nmf_metaprograms$clusters) %>% cumsum()
limits <- c(1, limits[-length(limits)])

x1 <- limits[seq(from = 1, to = length(limits) - 1, by = 1)]
x2 <- limits[seq(from = 2, to = length(limits), by = 1)]
x1 <- factor(levels(malignant_nmf_metaprograms$nmf_intersect_meltI$Var1)[x1],
             levels(malignant_nmf_metaprograms$nmf_intersect_meltI$Var1))
x2 <- factor(levels(malignant_nmf_metaprograms$nmf_intersect_meltI$Var1)[x2],
             levels(malignant_nmf_metaprograms$nmf_intersect_meltI$Var1))
y1 <- x1; y2 <- x2
f <- data.frame(x1 = x1, x2 = x2, y1 = y1, y2 = y2)

ggplot(data = malignant_nmf_metaprograms$nmf_intersect_meltI,
       aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
  geom_rect(data = f, aes(x = NULL,y = NULL, xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            color = "black", linetype = "dashed", fill = NA, size = .75, inherit.aes = F) +
  scale_y_discrete(limits = rev)
