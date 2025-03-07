## Glioblastoma Cellular Analysis of Resistance and Evolution consortium paper analyses

### Overview
The Cellular Analysis of Resistance and Evolution Glioblastoma (CARE GBM) consortium used both longitudinally collected single nucleus RNA sequencing and bulk DNA sequencing to redefine GBM (IDH-wild-type) cellular states, characterize intertumoral heterogeneity, and study tumor evolution. 

The code in this respository reflects analyses and figures associated with two separate manuscripts, "The Multi-Layered Transcriptional Architecture of Glioblastoma Ecosystems" (i.e., Nomura, Spitzer, Johnson, Garofano, et al) and "Deciphering the Longitudinal Trajectories of Glioblastoma by Integrative Single-Cell Genomics" (i.e., Spitzer, Johnson, Nomura, Garofano, et al). 

The `R` scripts directory reflects major scripts used to analyze the data associated with these manuscripts. R script file names are separated by paper using the listed first author's last name and figure panel to which the code pertains (e.g., `nomura_fig1a.R`).

### Data download
The CARE GBM processed snRNA patient-level 10x count matrices can be downloaded from [GSE274546](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546) and the processed cohort-level Smart-Seq2 gene expression matrix can be accessed at [GSE274548](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274548).

The CARE GBM DNA data can be downloaded from the `Tables` page [here](https://www.synapse.org/Synapse:syn26464346/tables/) and the `Files` page [here](https://www.synapse.org/Synapse:syn26464346/files/). It is also possible to query the data directly using the API by using queries. You can read more about that [here](https://docs.synapse.org/articles/tables.html).
