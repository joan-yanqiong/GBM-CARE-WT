#################################
## Title: Extended Data Figure 4b in Nomura et al
## Date: 2025.02.15
## Author: Kevin Johnson
## Description: Produce heatmap for signature scores correlation across published "functional signatures" (e.g., connectivity) and CARE 10x gene expression metaprogramss
#################################

# Necessary packages
library(tidyverse) # v1.3.1
# Custom minimalist plotting theme
source("R/plot_theme.R")

# Signature scores calculated (score_within_samples) for metaprograms and public signatures from Ivy GAP, Hai et al (PMID: 38320988), Taylor et al (PMID: 37914930), and Venkataramani et al. (PMID: 35914528).
mp_public_scores <- readRDS("data/nomura_edf4b_malignant_mp_public_scores.RDS")

mp_scores_mat <- mp_public_scores %>% 
  dplyr::select(MP_1_RP:IVY_LE) %>% 
  as.matrix()
cormat <- round(cor(mp_scores_mat), 2)

cormat_sig_df <- data.frame(cormat[1:15, 16:23])
cormat_sig_df$sig <- rownames(cormat_sig_df)
cormat_long <- cormat_sig_df %>% 
  pivot_longer(cols= c(TM_connected_Greenwald:IVY_LE),
               names_to = "Sig",
               values_to = "Cor") %>% 
  mutate(sig = recode(sig, `MP_1_RP` = "RP",
                      `MP_2_OPC` = "OPC-like",
                      `MP_3_CC` = "Cell cycle",
                      `MP_4_AC` = "AC-like",
                      `MP_5_Hypoxia` = "Hypoxia",
                      `MP_6_MES` = "MES-like",
                      `MP_7_NPC` = "NPC-like",
                      `MP_8_GPC` = "GPC-like",
                      `MP_9_ExN` = "NEU-like",
                      `MP_10_Stress1` = "Stress",
                      `MP_11_MIC` = "MIC",
                      `MP_12_LQ` = "LQ",
                      `MP_13_Cilia` = "Cilia",
                      `MP_14_NRGN` = "NRGN",
                      `MP_15_Stress2` = "Stress2")) %>% 
  mutate(Sig = recode(Sig, `Synaptic_Taylor` = "Synaptic (Taylor)",
                      `TM_connected_Greenwald` = "Connectivity (Hai)",
                      `TM_Taylor` = "Microtubes (Taylor)",
                      `invasive_Greenwald` = "Invasive (Venkataramani)"))
# Additional analysis.
# Produce plot for Ivy GAP signatures to show that leading edfe is enriched for scores in NEU-like, NPC-like metaprograms, cellular tumor is enriched for OPC/GPC-like. Pseudopallisading around necrosis is enriched for stress, MES-like, and Hypoxia.
cormat_long_order <- cormat_long %>% 
  filter(Sig=="IVY_LE") %>% 
  dplyr::arrange(Cor)

cormat_long <- cormat_long %>%
  mutate(sig = factor(sig, levels = rev(cormat_long_order$sig)),
         Sig = factor(Sig, levels = rev(c("IVY_LE", "IVY_CT", "IVY_CTmvp", "IVY_CTpan", "Connectivity (Hai)", "Microtubes (Taylor)","Invasive (Venkataramani)", "Synaptic (Taylor)"))))

ggplot(data = cormat_long %>% 
         filter(!Sig%in%c("Connectivity (Hai)", "Microtubes (Taylor)","Invasive (Venkataramani)", "Synaptic (Taylor)")), aes(y=Sig,x=sig, fill = Cor)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "dodgerblue", high = "#FF0000", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Signature score\ncorrelation") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1), 
        text=element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.position="top")+
  labs(y= "Public gene signatures", x= "Malignant metaprograms") 

# Sorting by correlation coefficient with connectivity to better visualize the difference in cell state scores across different functional categories.
cormat_long_order <- cormat_long %>% 
  filter(Sig== "Connectivity (Hai)") %>% 
  dplyr::arrange(Cor)

cormat_long <- cormat_long %>%
  mutate(sig = factor(sig, levels = rev(cormat_long_order$sig)),
         Sig = factor(Sig, levels = rev(c("IVY_LE", "IVY_CT", "IVY_CTmvp", "IVY_CTpan", "Connectivity (Hai)", "Microtubes (Taylor)","Invasive (Venkataramani)", "Synaptic (Taylor)")))) %>% 
  # Remove low-quality (e.g., RP and LQ metaprograms) and repetitive metaprograms (e.g., NRGN and Stress2).
  filter(!sig%in%c("RP", "NRGN", "Stress2", "LQ", "MIC"))

pdf(file = "figures/nomura_edf4b_malignant_connectivity_invasion_scores.pdf", height = 5, width = 7, bg = "transparent", useDingbats = FALSE)
ggplot(data = cormat_long %>% 
         filter(Sig%in%c("Connectivity (Hai)", "Microtubes (Taylor)","Invasive (Venkataramani)", "Synaptic (Taylor)")),
       aes(y=Sig,x=sig, fill = Cor)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "dodgerblue", high = "#FF0000", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Signature score\ncorrelation") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1), 
        text=element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.position="top")+
  labs(y= "Functional signatures", x= "GBM states") 
dev.off()


### END ###