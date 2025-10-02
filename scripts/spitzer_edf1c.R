#################################
## Title: Spitzer et al. Extended Data Figure 1c - genetic landscape across time points
## Date: 2025.02.18
## Author: Kevin Johnson
## Description: Generate a stacked bar plot for all mutations plus patient-specific visualization of altered genes
#################################

# Necessary packages
library(tidyverse) # v1.3.1
library(egg) # v0.4.5
library(grid) # base package
source("R/plot_theme.R")

# Enumerating mutations where there is coverage of at least 15x to determine relative longitudinal mutation fraction.
mut_freq_prop_case <- readRDS("data/genetic_patient_mutation_fraction.RDS")

# To determine arm gain/loss requires that a CNA segment be consistently gained or lost. Samples with disrupted segments where some proportions are gain while others are lost are assigned an NA.
arm_heatmap <- readRDS("data/genetic_patient_chr710_arm_status.RDS")
# Selected gene-level copy number alterations
cnv_heatmap <- readRDS("data/genetic_patient_gene_cna.RDS")
# Selected gene-level driver gene mutations
mut_heatmap <- readRDS("data/genetic_patient_gene_mutation.RDS")

# We'll be sorting this panel based on the fraction of shared mutations across time points.
sort_df <- mut_freq_prop_case %>%
  filter(mutation_type=="S") %>% 
  arrange(desc(mutation_percent))

patient_order <- unique(sort_df$patient)
snv_gene_order <- c("TERT", "PTEN", "EGFR", "TP53", "NF1", "RB1", "PIK3R1")
cnv_gene_order <- c("CDKN2A", "EGFR", "CDK4", "PDGFRA", "MDM2")
arm_order <- rev(c("7p", "7q", "10p", "10q"))

arm_heatmap <- arm_heatmap %>% mutate(patient = factor(patient, levels = patient_order),
                                      segment_qual = ifelse(!is.na(cnv_state), "Contiguous", "Disrupted"),
                                      cnv_state = ifelse(cnv_state == "neut", NA, cnv_state),
                                      arm = factor(arm, levels = arm_order))
cnv_heatmap <- cnv_heatmap %>% mutate(patient = factor(patient, levels = patient_order),
                                      gene_symbol = factor(gene_symbol, levels = cnv_gene_order))
mut_heatmap <- mut_heatmap %>% mutate(patient = factor(patient, levels = patient_order),
                                      gene_symbol = factor(gene_symbol, levels = snv_gene_order),
                                      covered = factor(covered, levels = c("Coverage < 5x", "Coverage >= 5x", "Coverage >= 15x", "Coverage >= 30x")))

### ### ### ### ### ### ###
## Common plotting elements
### ### ### ### ### ### ###

plot_grid     <- facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")
plot_theme    <- theme_bw(base_size = 14) + theme(axis.title = element_text(size = 14),
                                                  axis.text = element_text(size=14),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank())
null_legend   <- theme(legend.position = 'none')
null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
null_y        <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())
bottom_x      <- theme(axis.text.x=element_blank())
null_facet    <- theme(strip.background = element_blank(),
                       strip.text.x = element_blank())
top_margin    <- theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) 
middle_margin <- theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin <- theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))


plot_grid     <- facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")
plot_theme    <- theme_bw(base_size = 14) + theme(axis.title = element_text(size = 14),
                                                  axis.text = element_text(size=14),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank())
null_legend   <- theme(legend.position = 'none')
null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
null_y        <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank())
bottom_x      <- theme(axis.text.x=element_blank())
null_facet    <- theme(strip.background = element_blank(),
                       strip.text.x = element_blank())
top_margin    <- theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) 
middle_margin <- theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin <- theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))


gg_cbind <- function(..., widths = NULL) {
  if(length(match.call()) - 2 != length(widths))
    message("Number of widths does not match number of columns")
  gg <- gtable_cbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$widths[panels] <- unit(widths, "null")
  return(gg)
}


gg_rbind <- function(..., heights = NULL, ncol = 2) {
  if(length(match.call()) - 3 != length(heights))
    message("Number of heights does not match number of rows")
  gg <- gtable_rbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$heights[panels] <- unit(rep(heights,each = ncol), "null")
  return(gg)
}


### ### ### ### ### ### ###
## Plot proportions shared/private per patient
### ### ### ### ### ### ###
gg_mut_freq_prop_case <- ggplot(mut_freq_prop_case, aes(x = patient, y = mutation_percent, fill=mutation_type)) +
  geom_bar(stat="identity") +
  labs(y = "Proportion\nof mutations") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F"), na.value ="white",
                    labels=c('T2', 'T1', 'T1&T2')) +
  plot_grid

### ### ### ### ### ### ###
## Plot proportions shared/private per patient for mutation proportions, driver SNVs, and CNAs
### ### ### ### ### ### ###
gg_mut_freq_prop_case <- ggplot(mut_freq_prop_case, aes(x = patient, y = mutation_percent, fill=mutation_type)) +
  geom_bar(stat="identity") +
  labs(y = "Proportion\nof mutations") +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F"), na.value ="white",
                    labels=c('T2', 'T1', 'T1&T2')) +
  plot_grid


gg_mut_heatmap <- ggplot(mut_heatmap %>% 
           filter(gene_symbol%in%c("TERT","TP53", "PTEN", "EGFR", "NF1", "RB1", "PIK3R1")), 
         aes(x=patient, y=gene_symbol)) +
  geom_tile(aes(fill = variant_call)) +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F"), na.value ="white",
                    labels=c("R" = 'T2', "P"= 'T1', "S" = 'T1&T2'))  +
  labs(y = "Gene-level SNV", fill = "SNV evolution", x="")


gg_cnv_heatmap <-
  ggplot(cnv_heatmap %>% 
           filter(gene_symbol%in%c("CDKN2A", "EGFR", "CDK4", "PDGFRA", "MDM2")), 
         aes(x=patient, y=gene_symbol)) +
  geom_tile(aes(fill = cnv_change)) +
  scale_fill_manual(values=c("R" = "#2FB3CA", "P" ="#CA2F66", "S"="#CA932F"), na.value ="white",
                    labels=c("R" = 'T2', "P"= 'T1', "S" = 'T1&T2')) +
  labs(y = "Gene-level CNA", fill = "CNA evolution", x = "Patient") 


### ### ### ### ### ### ###
### Plot figure 1 heatmap
### ### ### ### ### ### ###
figms <- gg_rbind(gtable_frame(ggplotGrob(gg_mut_freq_prop_case + plot_grid + plot_theme + null_legend + null_x + null_facet)),
                  gtable_frame(ggplotGrob(gg_mut_heatmap + plot_grid + plot_theme + null_x + null_facet)),
                  gtable_frame(ggplotGrob(gg_cnv_heatmap + plot_grid + plot_theme + null_facet + bottom_margin + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))),
                  heights = c(4,5,4),
                  ncol = 1)
plot(figms)


pdf(file = "figures/spitzer_edf1c_longitudinal_genetics.pdf", height = 8, width = 10, bg = "transparent", useDingbats = FALSE)
grid.newpage()
grid.draw(figms)
dev.off()
