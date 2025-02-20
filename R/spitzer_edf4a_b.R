#################################
## Title: Wang et al sankey plots (Spitzer et al. EDF4a-b panels)
## Date: 2025.02.10
## Author: Kevin Johnson
## Description: Produce Sankey plots that characterize compositional and malignant group changes across time points.
#################################

# Necessary packages
library(networkD3) # v0.4
library(dplyr) # v1.1.4
library(htmlwidgets) # v1.6.4

# Diaz composition groups input data
links_comp <-read.delim("data/diaz_composition_groups_sankey_input.txt", sep = '\t', header = T)
links_comp$conservation <- ifelse(links_comp$source==links_comp$target, "conserved", "divergent")
table(links_comp$conservation)

# Retrieve the total amount of cases for whether the composition group remains the same ("conserved") or changes ("divergent")
links_comp %>% 
  group_by(conservation) %>% 
  summarise(total = sum(value))

# Total number of cases that shift compositional groups:
17/(17+8)

# Insert T1 and T2 nomenclature
links_comp$source <- paste0(links_comp$source, "_T1")
links_comp$target <- paste0(links_comp$target, "_T2")

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes_comp <- data.frame(
  name=c(as.character(links_comp$source), 
         as.character(links_comp$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe. We need to reformat it.
links_comp$IDsource <- match(links_comp$source, nodes_comp$name)-1 
links_comp$IDtarget <- match(links_comp$target, nodes_comp$name)-1

# Supply color code
my_color <- 'd3.scaleOrdinal() .domain(["HMF_T1", "IMF_T1","LMF-TAM_T1","LMF-OC_T1", "LMF-Mixed_T1", "HMF_T2", "IMF_T2","LMF-TAM_T2","LMF-OC_T2", "LMF-Mixed_T2"]) .range(["#FB8072", "#FDBFB8" , "#81B2D3","#B3DE69", "grey", "#FB8072", "#FDBFB8" , "#81B2D3","#B3DE69", "grey"])'

# Sankey plot for compositional clusters
p_comp <- sankeyNetwork(Links = links_comp, Nodes = nodes_comp, Source = "IDsource", Target = "IDtarget", 
                   Value = "value", NodeID = "name", colourScale=my_color,fontSize = 0,iterations = 0,nodeWidth = 90)


## Wang et al malignant group data
links_malignant <-read.delim("data/diaz_malignant_groups_sankey_input.txt", sep = '\t', header = T)
links_malignant$conservation <- ifelse(links_malignant$source==links_malignant$target, "conserved", "divergent")
table(links_malignant$conservation)

links_malignant %>% 
  group_by(conservation) %>% 
  summarise(total = sum(value))

# Total number of cases that have divergent malignant clusters across time points:
18/(18+7)

# Reformatting names for time points
links_malignant$source <- paste0(links_malignant$source, "_T1")
links_malignant$target <- paste0(links_malignant$target, "_T2")

# Repeat as above.
nodes_malignant <- data.frame(
  name=c(as.character(links_malignant$source), 
         as.character(links_malignant$target)) %>% unique()
)


links_malignant$IDsource <- match(links_malignant$source, nodes_malignant$name)-1 
links_malignant$IDtarget <- match(links_malignant$target, nodes_malignant$name)-1

# Difference in name and color structure for Sankey plots
my_color <- 'd3.scaleOrdinal() .domain(["AC_T1", "GPC_T1","MES/Hypoxia_T1","Mixed_T1", "OPC/NPC_T1", "AC_T2", "GPC_T2","MES/Hypoxia_T2","Mixed_T2", "OPC/NPC_T2"]) .range(["#66CC99", "#CC6633" , "#9999CC","#CC6699", "#FFCC33", "#66CC99", "#CC6633" , "#9999CC","#CC6699", "#FFCC33"])'

# Make the interactive network
p_malignant <- sankeyNetwork(Links = links_malignant, Nodes = nodes_malignant, Source = "IDsource", Target = "IDtarget", 
                   Value = "value", NodeID = "name", colourScale=my_color,fontSize = 0,iterations = 0,nodeWidth = 90)
p_malignant


# Save outputs as html for inspection.
# Export PDFs using Safari browser after resizing.
saveWidget(p_comp, file="figures/wang_longitudinal_compositional_groups.html")
saveWidget(p_malignant, file="figures/wang_longitudinal_malignant_groups.html")