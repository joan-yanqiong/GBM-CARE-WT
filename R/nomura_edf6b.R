#################################
## Title: Extended Data Figure 6 panel b in Nomura et al - Pathway enrichment analysis of PCA signatures
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: Produce a pathway enrichment graph of PCA signatures
#################################


ego_mps <- lapply(pca_gene_sigs, function(x) {
  clusterProfiler::enrichGO(gene          = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, x, 'ENTREZID', 'SYMBOL'),
                            universe      = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, names(rm), 'ENTREZID', 'SYMBOL'),
                            OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            minGSSize     = 50,
                            maxGSSize     = 300,
                            pool          = T,
                            readable      = TRUE)
})

####################################################################################################################################
# Plot panels
####################################################################################################################################

edo_mp <- enrichplot::pairwise_termsim(ego_mps[[1]])
enrichplot::emapplot(edo_mp, showCategory = 10)

edo_mp <- enrichplot::pairwise_termsim(ego_mps[[3]])
enrichplot::emapplot(edo_mp, showCategory = 10)

edo_mp <- enrichplot::pairwise_termsim(ego_mps[[5]])
enrichplot::emapplot(edo_mp, showCategory = 10)
