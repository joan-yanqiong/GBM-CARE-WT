library(dplyr)
library(tibble)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scales)
library(reshape2)
library(ggpubr)
library(ggridges)
library(infercna)
library(scalop)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
library(purrr)
library(DoubletFinder)
library(survival)
library(survminer)

source("R/GBM-CARE-WT_analysis_utils.R")
source("R/GBM-CARE-WT_CNA_utils.R")
source("R/GBM-CARE-WT_NMF.R")
source("R/GBM-CARE-WT_State_utils.R")

ANALYSIS_ROOT <- "~/"
DATA_ROOT <- paste0(ANALYSIS_ROOT, "data/")
TABLES_ROOT <- paste0(ANALYSIS_ROOT, "tables/")
SDS_SAMPLE_PATH <- paste0(DATA_ROOT, "sample_sds/")
NMSDS_SAMPLE_PATH <- paste0(DATA_ROOT, "sample_nmsds/")

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Load the per-sample meta-data
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

sample_data <- as_tibble(read.csv("path_to_meta_data_file/sample_data.csv", stringsAsFactors = F))

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Merge all individual count matrices to a per-sample list
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

path <- "path_to_count_matrices_location"

dirs_list <- list.files(path)

umi_data_all <- lapply(dirs_list, function(s) {
  message(s)
  umi_data <- Read10X(data.dir = paste0(path, s, "/outs/filtered_feature_bc_matrix/"))
  colnames(umi_data) <- paste0(s, "-", colnames(umi_data))
  umi_data
})
gc()

names(umi_data_all) <- dirs_list

saveRDS(object = umi_data_all, file = paste0(DATA_ROOT, "umi_data_list.RDS"))

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Prepare sample-data and meta-data files for analysis
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

sample_data$Patient_factor <- factor(sample_data$Patient, unique(sample_data$Patient))
sample_data$Timepoint_factor <- factor(sample_data$Timepoint, unique(sample_data$Timepoint))

sample_data$ID <- paste0(sample_data$Patient, sample_data$Timepoint)
sample_data$ID_extended <- paste0(sample_data$ID, "_", sample_data$Sample)

sample_data$ID_factor <- factor(sample_data$ID, unique(sample_data$ID))
sample_data$ID_extended_factor <- factor(sample_data$ID_extended, unique(sample_data$ID_extended))
table(sample_data$ID_factor)

colnames(sample_data)[grep("Processed", colnames(sample_data))] <- "Lab"

saveRDS(object = sample_data, file = paste0(DATA_ROOT, "sample_data.RDS"))

cellids <- unname(unlist(lapply(umi_data_all, colnames), recursive = F))
length(cellids)

meta_data <- tibble(CellID = cellids,
                    Sample = sapply(strsplit(CellID, split = "-"), function(y) y[[1]]))
table(meta_data$Sample)

meta_data <- meta_data %>% filter(Sample %in% sample_data$Sample)

meta_data <- meta_data %>% left_join(sample_data %>% select(Sample, Patient, Timepoint, PvsR, Patient_factor, Timepoint_factor, ID, ID_extended, ID_factor, ID_extended_factor,
                                                            Origin, Lab, Pathology, IDHstatus, MGMTstatus, SurgicalIntervalmonth, OSmonth, PatientStatus),
                                     by = "Sample")

meta_data <- as.data.frame(meta_data)
rownames(meta_data) <- meta_data$CellID

dim(meta_data)

saveRDS(object = meta_data, file = paste0(DATA_ROOT, "meta_data.RDS"))

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Basic QC
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

registered()

#############################################################################################################################################################
# Compute cell QC metrics (complexity, MT expression)
#############################################################################################################################################################

complexity <- plapply(X = unique(meta_data$Sample),
                      FUN = function(sname) {
                        m <- umi_data_all[[sname]]
                        m <- m != 0
                        colSums(m)
                      },
                      pack = T)

complexity <- complexity[meta_data$CellID]
length(complexity)

meta_data$Complexity <- complexity

meta_data %>%
  group_by(ID_factor, Sample, Origin, Lab) %>%
  summarise(N = n(),
            Mean = mean(Complexity),
            SD = sd(Complexity),
            Median = median(Complexity),
            Q25 = quantile(Complexity, .25),
            Q75 = quantile(Complexity, .75)) %>%
  View()

percent.mt <- plapply(X = unique(meta_data$Sample),
                      FUN = function(sname) {
                        m <- umi_data_all[[sname]]
                        mt_genes <- rownames(m)[grep("^MT-", rownames(m))]
                        mt_sum <- colSums(m[mt_genes, ])
                        all_sum <- colSums(m)
                        setNames((mt_sum / all_sum) * 100, colnames(m))
                      },
                      pack = T)

percent.mt <- percent.mt[meta_data$CellID]
length(percent.mt)

meta_data$Percent.MT <- percent.mt

meta_data %>%
  filter(Complexity >= 1000) %>%
  group_by(Patient_factor, Timepoint, Sample, Origin, Lab) %>%
  summarise(N = n(),
            Mean = mean(Percent.MT),
            SD = sd(Percent.MT),
            Median = median(Percent.MT),
            Q25 = quantile(Percent.MT, .25),
            Q75 = quantile(Percent.MT, .75),
            Q95 = quantile(Percent.MT, .95)) %>%
  View()

qc_stats <- meta_data %>%
  group_by(ID, Sample) %>%
  summarise(Ncells_PreQC = n(),
            MeanComplexity = mean(Complexity),
            MedianComplexity = median(Complexity),
            Q005Complexity = quantile(Complexity, .005),
            Q025Complexity = quantile(Complexity, .025),
            Q25Complexity = quantile(Complexity, .25),
            Q75Complexity = quantile(Complexity, .75),
            Q975Complexity = quantile(Complexity, .975))

meta_data <- meta_data %>% filter(Complexity >= 200 & Complexity <= 10000, Percent.MT <= 3)
dim(meta_data)

saveRDS(object = meta_data, file = paste0(DATA_ROOT, "meta_data.RDS"))

qc_stats <- left_join(x = meta_data %>% group_by(ID, Sample) %>% summarise(Ncells_PostQC = n()),
                      y = qc_stats,
                      by = c("ID", "Sample"))
qc_stats$PercentCellsDropped <- round((1 - qc_stats$Ncells_PostQC / qc_stats$Ncells_PreQC) * 100, 1)

qc_stats <- qc_stats %>% select(ID, Sample,
                                Ncells_PreQC, Ncells_PostQC, PercentCellsDropped,
                                MeanComplexity, MedianComplexity,
                                Q025Complexity, Q25Complexity, Q75Complexity, Q975Complexity)

qc_stats %>%
  arrange(ID) %>%
  write.table(file = paste0(TABLES_ROOT, "qc_stats_tbl.tsv"), sep = '\t', row.names = F, col.names = T)

#############################################################################################################################################################
# Select genes
#############################################################################################################################################################

selected_genes <- select_genes(umi_data_list = umi_data_all, md = meta_data, verbose = T, plot = T)

saveRDS(object = selected_genes, file = paste0(DATA_ROOT, "selected_genes.RDS"))

#############################################################################################################################################################
# Generate UMI matrix for seurat dataset
#############################################################################################################################################################

umi_data <- lapply(unique(meta_data$Sample), function(sname) {
  m <- umi_data_all[[sname]]
  m <- m[selected_genes, meta_data$CellID[meta_data$Sample == sname]]
  m
})
umi_data <- do.call(cbind, umi_data)
gc()
dim(umi_data)

saveRDS(object = umi_data, file = paste0(DATA_ROOT, "umi_data.RDS"))

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Create a top-level Seurat object
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

sds <- new_seurat_object(m = umi_data[, meta_data$CellID], md = meta_data, project = "GBM-CARE_integrated_analysis", verbose = T)
dim(sds)

DimPlot(sds, reduction = "umap", group.by = "Patient") +
  DimPlot(sds, reduction = "umap", group.by = "Timepoint") +
  DimPlot(sds, reduction = "umap", group.by = "Origin") +
  DimPlot(sds, reduction = "umap", group.by = "Lab")

DimPlot(sds, reduction = "umap", group.by = "IDHstatus") +
  DimPlot(sds, reduction = "umap", group.by = "MGMTstatus")

saveRDS(object = sds, file = paste0(DATA_ROOT, "sds.RDS"))

umap_res <- Reductions(sds, "umap")
meta_data$UMAP1 <- umap_res@cell.embeddings[meta_data$CellID, 1]
meta_data$UMAP2 <- umap_res@cell.embeddings[meta_data$CellID, 2]

saveRDS(object = meta_data, file = paste0(DATA_ROOT, "meta_data.RDS"))

rm(umap_res)
gc()

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Create the sample data object for sample specific cell type annotations, CNA inference and doublet detection
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

npcs <- 100

samples <- unique(meta_data$Sample)
length(samples)

start_time <- Sys.time()
for(i in 1:length(samples)) {
  
  gc()
  
  sname <- samples[i]
  
  message(paste0("******************************* ", sname, " - start, i = ", i, " *******************************"))
  
  cellids <- rownames(meta_data)[meta_data$Sample == sname]
  
  tmp <- new_seurat_object(m = umi_data[, cellids], md = meta_data[cellids, ],
                           project = paste0("GBM-Longitudinal_", sname), nfeatures = 3000, verbose = T)
  dim(tmp)
  
  
  print("Saving")
  st <- Sys.time()
  saveRDS(object = tmp, file = paste0(SDS_SAMPLE_PATH, "sds_", sname, ".RDS"))
  et <- Sys.time()
  print(et - st)
  
  message(paste0("******************************* ", sname, " - end, i = ", i, " *******************************"))
}
end_time <- Sys.time()
end_time - start_time

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Annotate cells for TME cell types
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

pathways <- load_c8_pathways()

macrophage <- derive_markers(pathways, "MACROPHAGE|MICROGLIA")
tcell <- derive_markers(pathways, "_T_CELL|_NK")
bcell <- derive_markers(pathways, "_B_CELL|_PLASMA", n_rep = 3)
endothel <- derive_markers(pathways, "ENDOTHEL")
pericyte <- derive_markers(pathways, "PERICYTE", n_rep = 2)
oligodendrocyte <- c("MOBP", "OPALIN", "CLDN11", "PLP1", "TF", "MBP", "MOG", "MAG", "EDIL3", "ST18", "CNTN2", "ANLN", "CNTN2",
                     "DOCK5", "TMEM144", "ENPP2", "SPOCK3", "PLEKHH1", "PIP4K2A", "BCAS1", "CERCAM",
                     "CNDP1", "CNTNAP4", "NKAIN2", "FRMD4B", "CLMN", "ZNF536", "AK5", "CTNNA3",
                     "ELMO1", "RNF220", "PEX5L", "SLC24A2", "PCSK6", "ST6GALNAC3", "IL1RAPL1", "MAP7",
                     "PLCL1", "UNC5C", "SLC5A11", "KLHL32", "UGT8", "CDK18", "ANO4", "FAM107B",
                     "GRM3", "MYRF", "TMEFF2", "FOLH1", "HHIP", "KCNH8", "PLD1")

sigs_list <- list(Macrophage = macrophage, Tcell = tcell, Bcell = bcell, Endothel = endothel, Pericyte = pericyte, Oligodendrocyte = oligodendrocyte)
lengths(sigs_list)

samples <- unique(meta_data$Sample)
length(samples)

# Annotate the cells
celltype_ann_step1 <- plapply(1:length(samples), function(i) {
  
  gc()
  
  sname <- samples[i]
  
  message(paste0("******************************* ", sname, " - start, i = ", i, " *******************************"))
  
  tmp <- load_sample(sname, SDS_SAMPLE_PATH)
  
  tmp <- AddModuleScore(object = tmp, features = sigs_list, name = "NormalSigs")
  
  d <- tmp@meta.data %>% select(CellID, Sample, starts_with("NormalSigs"))
  colnames(d)[grep("NormalSigs", colnames(d))] <- names(sigs_list)
  
  dm <- reshape2::melt(data = d,
                       id.vars = c("CellID", "Sample"))
  
  dm <- dm %>% group_by(CellID) %>% arrange(desc(value))
  dm <- dm[!duplicated(dm$CellID), ]
  dim(dm)
  dim(tmp)
  
  dm <- dm[dm$value > 1, ]
  dim(dm)
  
  new_idents <- setNames(rep("Presumed malignant", ncol(tmp)), colnames(tmp))
  new_idents[dm$CellID] <- as.character(dm$variable)
  
  Idents(tmp) <- new_idents
  
  colnames(tmp@meta.data[grep("NormalSigs", colnames(tmp@meta.data))]) <- names(sigs_list)
  
  d$CellType <- new_idents[d$CellID]
  
  message(paste0("******************************* ", sname, " - end, i = ", i, " *******************************"))
  
  return(d)
}, pack = F, profile = T, n_workers = 8)

celltype_ann_step1 <- do.call(rbind, celltype_ann_step1)

saveRDS(object = celltype_ann_step1, file = paste0(DATA_ROOT, "celltype_ann_step1.RDS"))

celltype_ann_step1$seurat_clusters <- sds$seurat_clusters[celltype_ann_step1$CellID]

celltype_ann_step1 %>%
  group_by(seurat_clusters) %>%
  summarise(Mean = mean(Oligodendrocyte), Q025 = quantile(Oligodendrocyte, .025), Q975 = quantile(Oligodendrocyte, .975)) %>%
  arrange(desc(Mean)) %>%
  View()

celltype_ann_step1_mlt <- melt(celltype_ann_step1)

ggplot(celltype_ann_step1_mlt, aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges(color = "black", scale = 1) +
  # scale_fill_viridis_c(name = "Score", option = "C", limits = c(0, .25), oob = squish) +
  xlab("Score") +
  ylab("CellType") +
  theme_gbm_pvsr(panel.border = element_blank())

celltype_ann_step1_mlt %>%
  group_by(variable, seurat_clusters) %>%
  summarise(Q25 = quantile(value, .25)) %>%
  group_by(variable) %>%
  arrange(desc(Q25)) %>%
  summarise(Q25 = first(Q25))

celltype_ann_step1_mlt <- celltype_ann_step1_mlt %>%
  filter(value > .4)
dim(celltype_ann_step1_mlt)

celltype_ann_step1_mlt <- celltype_ann_step1_mlt %>%
  group_by(CellID) %>%
  arrange(desc(value)) %>%
  summarise(CellID = first(CellID), CellType = first(variable))
table(celltype_ann_step1_mlt$CellType)

celltype_vec <- setNames(rep("Presumed malignant", ncol(sds)), colnames(sds))
celltype_vec[celltype_ann_step1_mlt$CellID] <- as.character(celltype_ann_step1_mlt$CellType)
table(celltype_vec)

Idents(sds) <- celltype_vec[colnames(sds)]

DimPlot(sds, reduction = "umap")

meta_data$CellType <- celltype_vec[meta_data$CellID]
table(meta_data$CellType)

ct_smp_tbl <- table(meta_data$Sample, meta_data$CellType) %>% as.matrix()
ct_smp_tbl <- ct_smp_tbl[, -6]
ct_smp_tbl[order(rowSums(ct_smp_tbl)), ]

rm(celltype_ann_step1)
rm(celltype_ann_step1_mlt)
rm(celltype_vec)
rm(ct_smp_tbl)
gc()

saveRDS(object = meta_data, file = paste0(DATA_ROOT, "meta_data.RDS"))

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# CNA classification algorithm 
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

#############################################################################################################################################################
# Compute CNA matrix per sample (defer correction step)
#############################################################################################################################################################

samples <- unique(meta_data$Sample)
length(samples)

cna_matrix_list <- lapply(samples, function(sname) {
  
  print(sname)
  
  md <- meta_data %>% filter(Sample == sname)
  
  m <- umi_data_all[[sname]]
  
  m <- umi2upm(m[, md$CellID])
  dim(m)
  
  infercna::useGenome("hg38")
  cna_matrix <- infercna::infercna(m = as.matrix(m), refCells = NULL, window = 100, isLog = F, n = 6000)
  dim(cna_matrix)
  cna_matrix <- t(cna_matrix)
  
  return(cna_matrix)
})

saveRDS(cna_matrix_list, paste0(DATA_ROOT, "cna_matrix_list.RDS"))

#############################################################################################################################################################
# Correct CNA scores using the reference populations
#############################################################################################################################################################

ref_types <- c("Oligodendrocyte", "Macrophage", "Endothel")
other_normal <- c("Tcell", "Bcell", "Pericyte")

cna_matrix_list_refcorrected <- lapply(cna_matrix_list, function(cna_matrix) {
  
  d <- meta_data %>%
    filter(CellID %in% rownames(cna_matrix))
  dim(d)
  
  print(paste0(unique(d$Sample)))
  
  ref_cells <- list(Oligodendrocyte = d$CellID[d$CellType_old == "Oligodendrocyte"],
                    Macrophage = d$CellID[d$CellType_old == "Macrophage"],
                    Endothel = d$CellID[d$CellType_old == "Endothel"])
  lengths(ref_cells)
  
  ref_cells <- ref_cells[lengths(ref_cells) >= 10]
  lengths(ref_cells)
  
  if(length(ref_cells) == 0) {
    print("No reference cells found - continuing to next one")
    return(NULL)  
  }
  
  print("Correcting CNA profiles using CNA values from <refCells>...")
  
  Args <- c(list(cna = t(cna_matrix), noise = 0.1, isLog = TRUE, zero_inflate = F), ref_cells)
  
  res <- do.call(refCorrect2, Args)
  
  return(t(res))
})
names(cna_matrix_list_refcorrected) <- samples

length(cna_matrix_list_refcorrected)
length(which(sapply(cna_matrix_list_refcorrected, is.null)))

rm(cna_matrix_list)
gc()

#############################################################################################################################################################
# TODO - compute CNA matrix for samples that don't have enough reference cells
#############################################################################################################################################################

null_samples <- names(cna_matrix_list_refcorrected)[sapply(cna_matrix_list_refcorrected, is.null)]

meta_data %>%
  filter(Sample %in% null_samples) %>%
  select(Sample, CellType_old) %>%
  table()

global_ref_cells <- meta_data %>%
  filter(IDHstatus == "WT", Sample %ni% null_samples, CellType_old %in% ref_types) %>%
  group_by(Sample, CellType_old) %>%
  filter(n() > 10) %>%
  sample_n(10)

saveRDS(global_ref_cells, paste0(DATA_ROOT, "global_ref_cells.RDS"))

table(global_ref_cells$Sample, global_ref_cells$CellType_old)
table(global_ref_cells$Sample, global_ref_cells$CellType_old) %>%
  rowSums() %>%
  sort()

m_ref <- lapply(names(umi_data_all), function(sname) {
  d <- global_ref_cells %>% filter(Sample == sname)
  if(nrow(d) == 0)
    return(NULL)
  m <- umi_data_all[[sname]][, d$CellID]
  return(m)
})
m_ref <- do.call(cbind, m_ref)
dim(m_ref)
m_ref[1:5, 1:5]

null_samples_cna_matrix_list <- lapply(null_samples, function(sname) {
  
  print(sname)
  
  md <- meta_data %>% filter(Sample == sname)
  
  md <- rbind(md, global_ref_cells)
  
  m <- umi_data_all[[sname]]
  
  m <- cbind(m, m_ref)
  
  m <- umi2upm(m[, md$CellID])
  dim(m)
  
  infercna::useGenome("hg38")
  cna_matrix <- infercna::infercna(m = as.matrix(m), refCells = NULL, window = 100, isLog = F, n = 6000)
  dim(cna_matrix)
  cna_matrix <- t(cna_matrix)
  
  return(cna_matrix)
})
names(null_samples_cna_matrix_list) <- null_samples

saveRDS(null_samples_cna_matrix_list, paste0(DATA_ROOT, "null_samples_cna_matrix_list.RDS"))

null_samples_cna_matrix_list_refcorrected <- lapply(names(null_samples_cna_matrix_list), function(sname) {
  
  print(sname)
  
  cna_matrix <- null_samples_cna_matrix_list[[sname]]
  
  d <- meta_data %>%
    filter(CellID %in% rownames(cna_matrix))
  dim(d)
  
  ref_cells <- list(Oligodendrocyte = d$CellID[d$CellType_old == "Oligodendrocyte"],
                    Macrophage = d$CellID[d$CellType_old == "Macrophage"],
                    # Tcell = d$CellID[d$CellType == "Tcell"],
                    Endothel = d$CellID[d$CellType_old == "Endothel"])
  lengths(ref_cells)
  
  ref_cells <- ref_cells[lengths(ref_cells) >= 10]
  lengths(ref_cells)
  
  if(length(ref_cells) == 0) {
    print("No reference cells found - continuing to next one")
    return(NULL)  
  }
  
  print("Correcting CNA profiles using CNA values from <refCells>...")
  
  Args <- c(list(cna = t(cna_matrix), noise = 0.1, isLog = TRUE, zero_inflate = F), ref_cells)
  
  res <- do.call(refCorrect2, Args)
  
  # Remove global reference cells
  sample_cells <- d %>% filter(Sample == sname)
  
  res <- res[, sample_cells$CellID]
  
  return(t(res))
})
names(null_samples_cna_matrix_list_refcorrected) <- null_samples

length(null_samples_cna_matrix_list_refcorrected)

lapply(null_samples_cna_matrix_list_refcorrected, dim)

# Merge null samples with all samples

cna_matrix_list_refcorrected[names(null_samples_cna_matrix_list_refcorrected)] <- null_samples_cna_matrix_list_refcorrected

length(cna_matrix_list_refcorrected)
length(which(sapply(cna_matrix_list_refcorrected, is.null)))

saveRDS(cna_matrix_list_refcorrected, paste0(DATA_ROOT, "cna_matrix_list_refcorrected.RDS"))

rm(cna_matrix_list)
rm(null_samples_cna_matrix_list)
gc()

#############################################################################################################################################################
# CNA classification STEP 0 - Pre-flight dataset setup
#############################################################################################################################################################

cna_sig_per_chr <- lapply(cna_matrix_list_refcorrected, function(cna_matrix) {
  
  if(is.null(cna_matrix))
    return(NULL)
  
  cna_genes <- splitGenes(colnames(cna_matrix))
  cna_genes <- cna_genes[-c(length(cna_genes))]
  names(cna_genes) <- paste0("CHR", names(cna_genes))
  
  cna_sig <- lapply(1:length(cna_genes), function(i) {
    tibble(CellID = rownames(cna_matrix),
           CHR = names(cna_genes)[i],
           Sig = rowMeans(cna_matrix[CellID, cna_genes[[i]]]))
  })
  cna_sig <- do.call(rbind, cna_sig)
  
  return(cna_sig)
})
cna_sig_per_chr <- do.call(rbind, cna_sig_per_chr)

cna_sig_per_chr$CHR <- factor(cna_sig_per_chr$CHR, levels = c(paste0("CHR", 1:22), "CHRX"))

cna_sig_per_chr <- cna_sig_per_chr %>%
  left_join(meta_data %>%
              select(CellID, Patient, Timepoint, Sample, Complexity, CellType_old, IDHstatus),
            by = "CellID")

# Limit the dataset to IDH-WT samples only
cna_sig_per_chr <- cna_sig_per_chr %>%
  filter(IDHstatus == "WT")
dim(cna_sig_per_chr)

unique(cna_sig_per_chr$Sample) %>% length()

setdiff(unique(meta_data$Sample[meta_data$IDHstatus == "WT"]), unique(cna_sig_per_chr$Sample))

cna_sig_per_chr <- cna_sig_per_chr %>%
  mutate(Group = case_when(CellType_old == "Presumed malignant" ~ "Presumed malignant",
                           CellType_old %in% ref_types ~ "Reference",
                           CellType_old %in% other_normal ~ "Other NM",
                           TRUE ~ "Unknown type")) %>%
  mutate(Group = factor(Group, c("Presumed malignant", "Reference", "Other NM")))

cna_sig_per_chr$isReference <- cna_sig_per_chr$Group == "Reference"
cna_sig_per_chr$isPN <- cna_sig_per_chr$Group == "Presumed malignant"

cna_sig_per_chr_stats <- cna_sig_per_chr %>%
  group_by(Sample, Group, CHR) %>%
  summarise(Mean = mean(Sig), Median = median(Sig), SD = sd(Sig), MAD = mad(Sig),
            Q01 = quantile(Sig, .01), Q99 = quantile(Sig, .99),
            Q25 = quantile(Sig, .25), Q75 = quantile(Sig, .75))

#############################################################################################################################################################
# CNA classification STEP 1 - Detect high frequency events
#############################################################################################################################################################

chr_sig_pm_vs_ref <- lapply(levels(cna_sig_per_chr_stats$CHR), function(chr) {
  x <- cna_sig_per_chr_stats %>% filter(CHR == chr, Group == "Presumed malignant") %>% pull(Mean)
  y <- cna_sig_per_chr_stats %>% filter(CHR == chr, Group == "Reference") %>% pull(Mean)
  tibble(CHR = chr, ES = median(x) - median(y), SignDiff = sign(median(x)) != sign(median(y)), pval = wilcox.test(x = x, y = y)$p.value)
})
chr_sig_pm_vs_ref <- do.call(rbind, chr_sig_pm_vs_ref)
chr_sig_pm_vs_ref$padj <- p.adjust(chr_sig_pm_vs_ref$pval, "holm")

ggboxplot(data = cna_sig_per_chr_stats, x = "CHR", y = "Mean", color = "Group", add = "jitter", add.params = list(size = .5),
          xlab = "", ylab = "Mean CNA signal per sample") +
  scale_color_manual(name = "", values = c("Presumed malignant" = "red", "Reference" = "dodgerblue", "Other NM" = "#00BA38")) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90)) +
  scale_y_continuous(limits = c(-.3, .3), oob = squish) +
  theme(panel.grid.major = element_line())

rm(cna_sig_per_chr_stats)
rm(chr_sig_pm_vs_ref)

#############################################################################################################################################################
# CNA classification STEP 1 - fit normal distribution to each chromosome
#############################################################################################################################################################

chr_sig_nd <- lapply(levels(cna_sig_per_chr$CHR), function(chr) {
  
  x <- cna_sig_per_chr %>%
    filter(isReference, CHR == chr) %>%
    pull(Sig)
  
  fit <- MASS::fitdistr(x, "normal")
  class(fit)
  
  para <- fit$estimate
  
  tibble(CHR = chr, Mean = para[1], SD = para[2])  
})
chr_sig_nd <- do.call(rbind, chr_sig_nd)

mean_vec <- setNames(chr_sig_nd$Mean, chr_sig_nd$CHR)
sd_vec <- setNames(chr_sig_nd$SD, chr_sig_nd$CHR)

norm_fit <- lapply(levels(cna_sig_per_chr$CHR), function(chr) tibble(CHR = chr,
                                                                     Sig = rnorm(n = 100000, mean = mean_vec[chr], sd = sd_vec[chr])))
norm_fit <- do.call(rbind, norm_fit)

d <- rbind(cna_sig_per_chr %>%
             select(CHR, Sig, Group),
           norm_fit %>%
             mutate(Group = "Classifier")) %>%
  mutate(Group = factor(Group, c("Presumed malignant", "Reference", "Other NM", "Classifier")))

ggplot(d %>% filter(CHR == "CHR10"), aes(x = Sig, y = after_stat(ncount), color = Group, linetype = Group)) +
  geom_freqpoly(binwidth = .01, show.legend = c(color = T, linetype = F), size = 1) +
  scale_color_manual(name = "", values = c("Reference" = "dodgerblue", "Presumed malignant" = "red",
                                           "Other NM" = "#00BA38", "Classifier" = "black")) +
  scale_linetype_manual(values = c("Reference" = "solid", "Presumed malignant" = "solid",
                                   "Other NM" = "solid", "Classifier" = "dashed")) +
  geom_vline(xintercept = mean_vec["CHR10"] - 1.96*sd_vec["CHR10"], linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = mean_vec["CHR10"] + 1.96*sd_vec["CHR10"], linetype = "dashed", color = "black", size = 1) +
  xlab("CHR10 signal") +
  ylab("Count (scaled to 1)") +
  scale_x_continuous(limits = c(-1, 1)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_gbm_pvsr(panel.border = element_blank())

chrs <- c("CHR1", "CHR4", "CHR6", "CHR7", "CHR8", "CHR10", "CHR14", "CHR21")

ggplot(d %>% filter(CHR %in% chrs), aes(x = Sig, y = after_stat(ncount), color = Group, linetype = Group)) +
  facet_wrap(~CHR, nrow = 2) +
  geom_freqpoly(binwidth = .01, show.legend = c(color = T, linetype = F), size = 1) +
  scale_color_manual(name = "", values = c("Reference" = "dodgerblue", "Presumed malignant" = "red",
                                           "Other NM" = "#00BA38", "Classifier" = "black")) +
  scale_linetype_manual(values = c("Reference" = "solid", "Presumed malignant" = "solid",
                                   "Other NM" = "solid", "Classifier" = "dashed")) +
  xlab("CNA signal") +
  ylab("Count (scaled to 1)") +
  scale_x_continuous(limits = c(-1, 1)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_gbm_pvsr(panel.spacing = unit(1, "lines"))

rm(d)

#############################################################################################################################################################
# CNA classification STEP 1 - Compute p-value per CHR for each cell
#############################################################################################################################################################

cna_sig_per_chr <- cna_sig_per_chr %>%
  mutate(CNAsig_z = (Sig - mean_vec[CHR]) / sd_vec[CHR])

# Two-sided test
cna_sig_per_chr <- cna_sig_per_chr %>%
  mutate(CNAsig_p = 2 * pnorm(abs(CNAsig_z), lower.tail = F))

ggplot(cna_sig_per_chr %>% filter(CHR == "CHR10"), aes(x = CNAsig_p, y = after_stat(count), fill = Group)) + 
  geom_histogram(binwidth = .01, position = "identity") +
  scale_fill_manual(name = "", values = c("Presumed malignant" = "red", "Reference" = "dodgerblue", "Other NM" = "#00BA38")) +
  geom_vline(xintercept = .05, linetype = "dashed", color = "black", size = 1) +
  xlab("P value") +
  ylab("Count") +
  # scale_y_log10() +
  theme_gbm_pvsr(panel.border = element_blank())

ggplot(cna_sig_per_chr %>% filter(CHR %in% chrs), aes(x = CNAsig_p, y = after_stat(count), fill = Group)) + 
  facet_wrap(~CHR, nrow = 2) +
  geom_histogram(binwidth = .01, position = "identity") +
  scale_fill_manual(name = "", values = c("Reference" = "dodgerblue", "Presumed malignant" = "red", "Other NM" = "#00BA38")) +
  geom_vline(xintercept = .05, linetype = "dashed", color = "black", size = 1) +
  xlab("P value") +
  ylab("Count") +
  scale_x_continuous(breaks = c(0, .5, 1)) +
  theme_gbm_pvsr(panel.spacing = unit(1, "lines")) +
  theme(panel.grid.major = element_line())

#############################################################################################################################################################
# CNA classification STEP 1 - FDR per CHR
#############################################################################################################################################################

cna_sig_per_chr <- cna_sig_per_chr %>%
  group_by(CHR) %>%
  mutate(CNAsig_padj = p.adjust(CNAsig_p, "fdr")) %>%
  ungroup()

cna_sig_per_chr$CNAsig_p_isSig <- cna_sig_per_chr$CNAsig_p < .05
cna_sig_per_chr %>%
  filter(CHR == "CHR10") %>%
  select(Group, CNAsig_p_isSig) %>%
  table()

cna_sig_per_chr$CNAsig_padj_isSig <- cna_sig_per_chr$CNAsig_padj < .05
tbl <- cna_sig_per_chr %>%
  filter(CHR == "CHR10") %>%
  select(Group, CNAsig_padj_isSig) %>%
  table()
tbl
tbl / rowSums(tbl) * 100

cna_sig_per_chr$CNAsig_FDR5 <- cna_sig_per_chr$CNAsig_padj < .05
cna_sig_per_chr$CNAsig_FDR1 <- cna_sig_per_chr$CNAsig_padj < .01

tbl <- cna_sig_per_chr %>%
  filter(CHR == "CHR10") %>%
  select(Group, CNAsig_FDR5) %>%
  table()
tbl / rowSums(tbl) * 100
tbl <- cna_sig_per_chr %>%
  filter(CHR == "CHR10") %>%
  select(Group, CNAsig_FDR1) %>%
  table()
tbl / rowSums(tbl) * 100

cna_events_per_chr <- cna_sig_per_chr %>%
  filter(CHR == "CHR10") %>%
  group_by(Group) %>%
  summarise(N = n(), Uncontrolled = sum(CNAsig_p_isSig) / N, FDR5 = sum(CNAsig_FDR5) / N, FDR1 = sum(CNAsig_FDR1) / N)

tbl <- melt(cna_events_per_chr %>% select(-N))

ggplot(tbl, aes(x = Group, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = paste0(round(value * 100, 1), "%"), hjust = "center", vjust = "center"),
            size = 7, color = "black", position = position_dodge(width = .75)) +
  scale_fill_brewer(name = "", palette = "Set1", labels = c("Uncontrolled", "5% FDR", "1% FDR")) +
  xlab("") +
  ylab("% cells classified with CHR10 event") +
  geom_hline(yintercept = .05, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = .01, linetype = "dashed", color = "black", size = 1) +
  scale_y_continuous(labels = percent)+
  theme_gbm_pvsr(panel.border = element_blank())

cna_sig_per_chr <- cna_sig_per_chr %>%
  mutate(CNA = case_when(Sig > 0 ~ "Amp",
                         Sig < 0 ~ "Del",
                         TRUE ~ "Balanced"))
table(cna_sig_per_chr$CHR, cna_sig_per_chr$CNA)

cna_events_per_chr <- cna_sig_per_chr %>%
  group_by(Group, CHR, CNA) %>%
  summarise(Uncontrolled = sum(CNAsig_p_isSig), FDR5 = sum(CNAsig_FDR5), FDR1 = sum(CNAsig_FDR1)) %>%
  left_join(cna_sig_per_chr %>%
              group_by(Group, CHR) %>%
              summarise(N = n()), by = c("Group", "CHR")) %>%
  mutate(Uncontrolled = Uncontrolled / N, FDR5 = FDR5 / N, FDR1 = FDR1 / N)

tbl <- melt(cna_events_per_chr %>% select(-N) %>% filter(CHR == "CHR10"))

ggplot(tbl, aes(x = Group, y = value, fill = variable)) +
  facet_wrap(~CNA) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = paste0(round(value * 100, 1), "%"), hjust = "center", vjust = "center"),
            size = 6, color = "black", position = position_dodge(width = .75)) +
  scale_fill_brewer(name = "", palette = "Set1", labels = c("Uncontrolled", "5% FDR", "1% FDR")) +
  xlab("") +
  ylab("% cells classified with CHR10 event") +
  geom_hline(yintercept = .05, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = .01, linetype = "dashed", color = "black", size = 1) +
  scale_y_continuous(labels = percent)+
  theme_gbm_pvsr()

cna_sig_per_chr <- cna_sig_per_chr %>%
  mutate(CNA = case_when(CNAsig_padj_isSig == T & Sig > 0 ~ 1,
                         CNAsig_padj_isSig == T & Sig < 0 ~ -1,
                         TRUE ~ 0))
table(cna_sig_per_chr$CHR, cna_sig_per_chr$CNA)

cna_sig_per_chr <- cna_sig_per_chr %>%
  mutate(CNAchr = case_when(CNA == 1 ~ "Amp",
                            CNA == -1 ~ "Del",
                            TRUE ~ "Balanced")) %>%
  mutate(CNAfct = factor(CNAchr, c("Del", "Balanced", "Amp")))

ggplot(cna_sig_per_chr %>% filter(CHR == "CHR10"), aes(x = CNAsig_z, y = after_stat(ncount), color = Group, linetype = Group)) +
  geom_freqpoly(binwidth = .1, show.legend = c(color = T, linetype = F), size = 1) +
  scale_color_manual(name = "", values = c("Reference" = "dodgerblue", "Presumed malignant" = "red", "Other NM" = "#00BA38")) +
  scale_linetype_manual(values = c("Reference" = "solid", "Presumed malignant" = "solid", "Other NM" = "solid")) +
  geom_vline(xintercept = zcrit_p_chr10, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = zcrit_fdr5_chr10, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = zcrit_fdr1_chr10, linetype = "dashed", color = "black", size = 1) +
  xlab("CHR10 signal [Z score]") +
  ylab("Count (scaled to 1)") +
  scale_x_continuous(limits = c(-5, 0)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_gbm_pvsr(panel.border = element_blank())

cna_sig_per_chr_stats <- cna_sig_per_chr %>%
  group_by(CHR, Group, CNA, CNAfct) %>%
  summarise(n = n()) %>%
  group_by(CHR, Group) %>%
  mutate(N = sum(n), Freq = n / N)
cna_sig_per_chr_stats$Freq <- cna_sig_per_chr_stats$Freq * sign(cna_sig_per_chr_stats$CNA)

ggplot(cna_sig_per_chr_stats %>% filter(CNA != "0"), aes(x = CHR, y = Freq, fill = CNAfct)) +
  facet_wrap(~Group) +
  geom_bar(color = "black", stat = "identity") +
  xlab("") +
  ylab("% cells with event") +
  scale_fill_manual(name = "CNA event", values = c("dodgerblue", "red")) +
  scale_y_continuous(labels = c("50%", "40%", "30%", "20%", "10%", "0%", "10%", "20%", "30%"), breaks = seq(-.5, .3, .1), limits = c(-.6, .4)) +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())

#############################################################################################################################################################
# CNA classification STEP 1 - Classify real events
#############################################################################################################################################################

WES_pairs_tbl <- read.csv(paste0(TABLES_ROOT, "analysis_tumor_pairs.csv"), header = T, stringsAsFactors = F)
WES_pairs_tbl <- WES_pairs_tbl[, -c(1:2)]
WES_pairs_tbl <- as_tibble(WES_pairs_tbl)

WES_pairs_tbl <- rbind(tibble(tumor_pair_barcode = WES_pairs_tbl$tumor_pair_barcode, case_barcode = WES_pairs_tbl$case_barcode,
                              tumor_barcode = WES_pairs_tbl$tumor_barcode_a, ID = WES_pairs_tbl$avishay_id_a),
                       tibble(tumor_pair_barcode = WES_pairs_tbl$tumor_pair_barcode, case_barcode = WES_pairs_tbl$case_barcode,
                              tumor_barcode = WES_pairs_tbl$tumor_barcode_b, ID = WES_pairs_tbl$avishay_id_b))
dim(WES_pairs_tbl)

d <- meta_data %>%
  filter(IDHstatus == "WT") %>%
  select(Patient, Sample, Timepoint, ID) %>%
  filter(!duplicated(ID)) %>%
  as_tibble()
dim(d)

WES_pairs_tbl <- WES_pairs_tbl %>%
  filter(ID %in% d$ID)
dim(WES_pairs_tbl)

WES_pairs_tbl <- WES_pairs_tbl %>%
  left_join(d, by = "ID")
dim(WES_pairs_tbl)

WES_CNA_calls <- read.csv(paste0(TABLES_ROOT, "analysis_gatk_cnv_by_arm.csv"), header = T, stringsAsFactors = F)
WES_CNA_calls <- WES_CNA_calls[, -c(1:2)]

WES_CNA_calls <- as_tibble(WES_CNA_calls)

WES_CNA_calls <- WES_CNA_calls %>%
  rename(tumor_barcode = aliquot_barcode) %>%
  filter(tumor_barcode %in% WES_pairs_tbl$tumor_barcode)
dim(WES_CNA_calls)

WES_CNA_calls <- WES_CNA_calls %>%
  left_join(WES_pairs_tbl %>%
              select(tumor_barcode, Patient, Timepoint, Sample, ID),
            by = "tumor_barcode")

colnames(WES_CNA_calls)[2] <- "CHR"

WES_CNA_calls$arm <- factor(WES_CNA_calls$arm, unique(WES_CNA_calls$arm))
WES_CNA_calls$CHR <- paste0("CHR", WES_CNA_calls$CHR)
WES_CNA_calls$CHR <- factor(WES_CNA_calls$CHR, unique(WES_CNA_calls$CHR))

WES_CNA_calls$arm_call[is.na(WES_CNA_calls$arm_call)] <- 0

WES_CNA_calls_chr <- WES_CNA_calls %>%
  group_by(tumor_barcode, Patient, Timepoint, Sample, CHR) %>%
  summarise(Call = sum(arm_call) / abs(sum(arm_call))) %>%
  ungroup() %>%
  mutate(Call = ifelse(is.na(Call), 0, Call))
table(WES_CNA_calls_chr$CHR, WES_CNA_calls_chr$Call)

WES_CNA_calls_chr$fDetected <- sapply(1:nrow(WES_CNA_calls_chr), function(i) {
  sample <- WES_CNA_calls_chr$Sample[i]
  chr <- WES_CNA_calls_chr$CHR[i]
  chr_call <- WES_CNA_calls_chr$Call[i]
  
  d <- cna_sig_per_chr %>%
    filter(Sample == sample, CHR == as.character(chr), isPN)
  
  if(nrow(d) == 0)
    return(0)
  
  res <- 0
  if(chr_call != 0) {
    res <- sum(d$CNA == chr_call & d$CNAsig_padj_isSig)
  } else {
    res <- sum(d$CNA == chr_call)
  }
  res <- res / nrow(d)
  
  return(res)
})

WES_CNA_calls_chr <- WES_CNA_calls_chr %>%
  mutate(Call_chr = case_when(Call == 0 ~ "Balanced",
                              Call == 1 ~ "Amp",
                              Call == -1 ~ "Del",
                              TRUE ~ "NA"))

ggboxplot(data = WES_CNA_calls_chr %>% filter(Call != 0), x = "CHR", y = "fDetected", add = "jitter",
          xlab = "", ylab = "% cells per sample with event") +
  facet_grid(rows = vars(Call_chr)) +
  scale_y_continuous(labels = percent) +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90), panel.spacing = unit(1, "lines")) +
  theme(panel.grid.major = element_line())

WES_CNA_calls_chr %>%
  group_by(CHR, Call, Call_chr) %>%
  summarise(n = n()) %>%
  group_by(CHR) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  filter(Call != 0) %>%
  ggline(x = "CHR", y = "Freq", color = "Call_chr", size = 1, xlab = "",
         ylab = "% samples with event\n(detected by WES)") +
  geom_hline(yintercept = .1, linetype = "dashed", color = "black", size = 1) +
  scale_color_brewer(name = "", palette = "Set1") +
  scale_y_continuous(labels = percent, breaks = c(0, .05, .1, .25, .5, .65)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90), panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

WES_CNA_calls_chr_stats <- WES_CNA_calls_chr %>%
  group_by(CHR, Call, Call_chr) %>%
  summarise(Mean = mean(fDetected))

WES_CNA_calls_chr_stats %>%
  filter(Call != 0) %>%
  ggline(x = "CHR", y = "Mean", color = "Call_chr", size = 1, xlab = "",
         ylab = "Average % cells with event\n(given detected by WES)") +
  scale_color_brewer(name = "", palette = "Set1") +
  scale_y_continuous(labels = percent, breaks = c(0, .05, .1, .2, .3, .4, .5)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90), panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

cna_sig_per_chr %>%
  group_by(Sample, CHR, CNA, CNAchr) %>%
  summarise(n = n()) %>%
  group_by(Sample, CHR) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  filter(CNA != 0) %>%
  group_by(CHR, CNA, CNAchr) %>%
  summarise(n = sum(Freq > .1), N = n(), Freq = n / N) %>%
  ggline(x = "CHR", y = "Freq", color = "CNAchr", size = 1, xlab = "",
         ylab = "% samples with event") +
  geom_hline(yintercept = .05, linetype = "dashed", color = "black", size = 1) +
  scale_color_brewer(name = "", palette = "Set1") +
  scale_y_continuous(labels = percent, breaks = c(0, .05, .1, .25, .5, .65)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90), panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

WES_CNA_calls_chr_balanced_discordant <- WES_CNA_calls_chr %>% filter(Call == 0)

WES_CNA_calls_chr_balanced_discordant$fDiscordant <- 1 - WES_CNA_calls_chr_balanced_discordant$fDetected

ggboxplot(data = WES_CNA_calls_chr_balanced_discordant, x = "CHR", y = "fDiscordant", add = "jitter",
          xlab = "", ylab = "% cells per sample with event") +
  geom_hline(yintercept = .05, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = .1, linetype = "dashed", color = "black", size = 1) +
  scale_y_continuous(labels = percent) +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90), panel.spacing = unit(1, "lines")) +
  theme(panel.grid.major = element_line())

wes_high_freq_cna_events <- WES_CNA_calls_chr %>%
  group_by(CHR, Call) %>%
  summarise(n = n()) %>%
  group_by(CHR) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  filter(Freq > .1, Call != 0)
wes_high_freq_cna_events$Event <- paste0(wes_high_freq_cna_events$CHR, "-", ifelse(wes_high_freq_cna_events$Call == 1, "Amp", "Del"))

snrnaseq_high_freq_cna_events <- cna_sig_per_chr %>%
  group_by(Sample, CHR, CNA) %>%
  summarise(n = n()) %>%
  group_by(Sample, CHR) %>%
  mutate(N = sum(n), Freq = n / N) %>%
  filter(CNA != 0) %>%
  group_by(CHR, CNA) %>%
  summarise(n = sum(Freq > .1), N = n(), Freq = n / N) %>%
  filter(Freq > .05)
snrnaseq_high_freq_cna_events$Event <- paste0(snrnaseq_high_freq_cna_events$CHR, "-", ifelse(snrnaseq_high_freq_cna_events$CNA == 1, "Amp", "Del"))

consensus_cna_events <- intersect(wes_high_freq_cna_events$Event, snrnaseq_high_freq_cna_events$Event)

wes_events_per_sample_tbl <- WES_CNA_calls_chr %>%
  filter(Call != 0)
wes_events_per_sample_tbl$Event <- paste0(wes_events_per_sample_tbl$CHR, "-", ifelse(wes_events_per_sample_tbl$Call == 1, "Amp", "Del"))

cna_sig_per_chr <- cna_sig_per_chr %>%
  mutate(Event = paste0(CHR, "-", CNAchr))

cna_sig_per_chr <- cna_sig_per_chr %>%
  left_join(WES_CNA_calls_chr %>%
              select(Sample, CHR, Call),
            by = c("Sample", "CHR"))
which(is.na(cna_sig_per_chr$Call)) %>% length()

no_call_samples <- setdiff(unique(cna_sig_per_chr$Sample), unique(wes_events_per_sample_tbl$Sample))

cna_sig_per_chr$Call[cna_sig_per_chr$Sample %in% no_call_samples] <- NA

cna_sig_per_chr$Call[is.na(cna_sig_per_chr$Call) & cna_sig_per_chr$CHR == "CHRX"] <- 0

cna_sig_per_chr$Call[is.na(cna_sig_per_chr$Call)] <- ifelse(cna_sig_per_chr$Event[is.na(cna_sig_per_chr$Call)] %in% consensus_cna_events, cna_sig_per_chr$CNA[is.na(cna_sig_per_chr$Call)], 0)

cna_sig_per_chr$isRealEvent <- cna_sig_per_chr$CNA & cna_sig_per_chr$Call
table(cna_sig_per_chr$CHR, cna_sig_per_chr$isRealEvent)
length(which(is.na(cna_sig_per_chr$isRealEvent)))

cna_sig_per_chr %>%
  group_by(Sample, CellID) %>%
  summarise(hasEvent = sum(isRealEvent) > 0) %>%
  group_by(Sample) %>%
  summarise(n = sum(hasEvent)) %>%
  filter(n == 0)

cna_sig_per_chr$CNAsign_chr <- ifelse(sign(cna_sig_per_chr$Sig) > 0, "Amp", "Del")

cna_sig_per_chr <- cna_sig_per_chr %>%
  mutate(Call_chr = case_when(Call == 1 ~ paste0(CHR, "-", "Amp"),
                              Call == -1 ~ paste0(CHR, "-", "Del"),
                              TRUE ~ paste0(CHR, "-", "Balanced")))
table(cna_sig_per_chr$CHR, cna_sig_per_chr$Call_chr)

cna_sig_per_chr <- cna_sig_per_chr %>%
  mutate(isSuspiciousEvent = !CNAsig_padj_isSig & CNAsig_p_isSig & (paste0(CHR, "-", CNAsign_chr) == Call_chr))
table(cna_sig_per_chr$CHR, cna_sig_per_chr$isSuspiciousEvent)

cna_sig_per_chr %>%
  group_by(CellID, CellType_old) %>%
  summarise(hasSuspicious = sum(isSuspiciousEvent) > 0) %>%
  ungroup() %>%
  select(CellType_old, hasSuspicious) %>%
  table()

rm(d)
rm(WES_CNA_calls_chr_balanced_discordant)
rm(WES_CNA_calls_chr_stats)
gc()

#############################################################################################################################################################
# CNA classification STEP 1 - Classify cells based on events
#############################################################################################################################################################

nm_cells <- cna_sig_per_chr %>%
  filter(Group %in% c("Reference", "Other NM")) %>%
  group_by(CellID, CellType_old) %>%
  summarise(nEvents = sum(isRealEvent), nSuspicipous = sum(isSuspiciousEvent)) %>%
  mutate(CNAclassification_step1 = ifelse(nEvents == 0 & nSuspicipous == 0, "Non-malignant", "Unresolved")) %>%
  ungroup()
table(nm_cells$CNAclassification_step1)
table(nm_cells$CellType_old, nm_cells$CNAclassification_step1)
table(nm_cells$CellType_old, nm_cells$CNAclassification_step1) / rowSums(table(nm_cells$CellType_old, nm_cells$CNAclassification_step1))
table(nm_cells$CNAclassification_step1, nm_cells$nEvents)
table(nm_cells$CNAclassification_step1, nm_cells$nSuspicipous)

pm_cells <- cna_sig_per_chr %>%
  filter(Group %ni% c("Reference", "Other NM")) %>%
  group_by(CellID, CellType_old, CHR) %>%
  mutate(isCHR7 = (CHR == "CHR7") & isRealEvent, isCHR10 = (CHR == "CHR10") & isRealEvent) %>%
  group_by(CellID, CellType_old) %>%
  summarise(nEvents = sum(isRealEvent), isCHR7 = sum(isCHR7) > 0, isCHR10 = sum(isCHR10) > 0) %>%
  mutate(CNAclassification_step1 = ifelse(nEvents == 0, "Non-malignant", "Malignant")) %>%
  ungroup()
table(pm_cells$CNAclassification_step1)
table(pm_cells$CNAclassification_step1, pm_cells$nEvents)
table(pm_cells$CNAclassification_step1, pm_cells$isCHR7)
table(pm_cells$CNAclassification_step1, pm_cells$isCHR10)

pm_cells %>%
  filter(CNAclassification_step1 == "Malignant") %>%
  select(isCHR7, nEvents) %>%
  table()

pm_cells %>%
  filter(CNAclassification_step1 == "Malignant") %>%
  select(isCHR10, nEvents) %>%
  table()

cna_sig_per_chr %>%
  filter(CHR == "CHR7", isPN) %>%
  select(isRealEvent) %>%
  table()

pm_nm_cells <- cna_sig_per_chr %>%
  filter(CellID %in% pm_cells$CellID[pm_cells$CNAclassification_step1 == "Non-malignant"])
dim(pm_nm_cells)
length(unique(pm_nm_cells$CellID))

pm_nm_cells_stats <- pm_nm_cells %>%
  group_by(CellID, CellType_old, CHR) %>%
  filter(!isRealEvent) %>%
  group_by(CellID, CellType_old) %>%
  summarise(nEvents = sum(CNAsig_padj_isSig))
table(pm_nm_cells_stats$nEvents)
table(pm_nm_cells_stats$nEvents) / nrow(pm_nm_cells_stats)

pm_nm_cells_stats <- pm_nm_cells %>%
  group_by(CellID) %>%
  summarise(nSuspicipous = sum(isSuspiciousEvent))
table(pm_nm_cells_stats$nSuspicipous)
table(pm_nm_cells_stats$nSuspicipous) / nrow(pm_nm_cells_stats)

pm_nm_cells <- pm_nm_cells %>%
  group_by(CellID, CellType_old) %>%
  summarise(nSuspicipous = sum(isSuspiciousEvent)) %>%
  mutate(CNAclassification_step1 = ifelse(nSuspicipous == 0, "Non-malignant", "Unresolved"))
table(pm_nm_cells$CNAclassification_step1)

one_event_cells <- pm_nm_cells %>%
  filter(nSuspicipous == 1) %>%
  pull(CellID)

one_event_cells <- cna_sig_per_chr %>%
  filter(CellID %in% one_event_cells, isSuspiciousEvent)
table(one_event_cells$CHR)
table(one_event_cells$CHR) / nrow(one_event_cells)

cell_classification_data <- rbind(nm_cells %>%
                                    select(CellID, CellType_old, CNAclassification_step1),
                                  pm_cells %>%
                                    filter(CNAclassification_step1 == "Malignant") %>%
                                    select(CellID, CellType_old, CNAclassification_step1),
                                  pm_nm_cells %>%
                                    select(CellID, CellType_old, CNAclassification_step1))
table(cell_classification_data$CellType_old, cell_classification_data$CNAclassification_step1)

cna_sig_per_chr <- cna_sig_per_chr %>%
  ungroup() %>%
  left_join(cell_classification_data %>%
              select(CellID, CNAclassification_step1),
            by = "CellID")

cna_sig_per_chr <- cna_sig_per_chr %>%
  mutate(GroupClass = case_when(isPN & CNAclassification_step1 == "Malignant" ~ "Malignant",
                                isPN & CNAclassification_step1 == "Non-malignant" ~ "PM-NM",
                                isPN & CNAclassification_step1 == "Unresolved" ~ "PM-U",
                                isReference & CNAclassification_step1 == "Non-malignant" ~ "REF-NM",
                                isReference & CNAclassification_step1 == "Unresolved" ~ "REF-U",
                                (Group == "Other NM") & CNAclassification_step1 == "Non-malignant" ~ "ONM-NM",
                                (Group == "Other NM") & CNAclassification_step1 == "Unresolved" ~ "ONM-U",
                                TRUE ~ "NA")) %>%
  mutate(GroupClass = factor(GroupClass, c("Malignant", "PM-NM", "REF-NM", "ONM-NM", "PM-U", "REF-U", "ONM-U")))
table(cna_sig_per_chr$GroupClass)

cna_sig_per_chr_stats <- cna_sig_per_chr %>%
  group_by(Sample, CHR, GroupClass) %>%
  summarise(Mean = mean(Sig))

ggboxplot(data = cna_sig_per_chr_stats, x = "CHR", y = "Mean", color = "GroupClass", add = "jitter", add.params = list(size = .5),
          xlab = "", ylab = "Mean CNA signal per sample") +
  scale_color_brewer(name = "", palette = "Set1") +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90)) +
  scale_y_continuous(limits = c(-.3, .3), oob = squish) +
  theme(panel.grid.major = element_line())

ggplot(cna_sig_per_chr %>% filter(CHR == "CHR10"), aes(x = Sig, y = after_stat(ncount), color = GroupClass)) +
  geom_freqpoly(binwidth = .01, show.legend = c(color = T, linetype = F), size = 1) +
  geom_freqpoly(data = norm_fit %>% filter(CHR == "CHR10"), mapping = aes(x = Sig, y = after_stat(ncount)),
                binwidth = .01, color = "black", linetype = "dashed", size = 1, inherit.aes = F) +
  scale_color_brewer(name = "Complexity threshold", palette = "Set1") +
  geom_vline(xintercept = mean_vec["CHR10"] - 1.96*sd_vec["CHR10"], linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = mean_vec["CHR10"] + 1.96*sd_vec["CHR10"], linetype = "dashed", color = "black", size = 1) +
  xlab("CHR10 signal") +
  ylab("Count (scaled to 1)") +
  scale_x_continuous(limits = c(-1, 1)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_gbm_pvsr(panel.border = element_blank())

ggplot(cna_sig_per_chr %>% filter(CHR %in% chrs), aes(x = Sig, y = after_stat(ncount), color = GroupClass)) +
  facet_wrap(~CHR, nrow = 2) +
  geom_freqpoly(binwidth = .01, show.legend = c(color = T, linetype = F), size = 1) +
  geom_freqpoly(data = norm_fit %>% filter(CHR %in% chrs), mapping = aes(x = Sig, y = after_stat(ncount)),
                binwidth = .01, color = "black", linetype = "dashed", size = 1, inherit.aes = F) +
  scale_color_brewer(name = "Complexity threshold", palette = "Set1") +
  xlab("CNA signal") +
  ylab("Count (scaled to 1)") +
  scale_x_continuous(limits = c(-1, 1)) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_gbm_pvsr(panel.border = element_blank())

cell_classification_data <- cell_classification_data %>%
  left_join(meta_data %>%
              select(CellID, CellType),
            by = "CellID")

cell_classification_data <- cell_classification_data %>%
  mutate(isPN = CellType_old == "Presumed malignant",
         isReference = CellType_old %in% ref_types,
         isONM = CellType_old %in% other_normal) %>%
  mutate(GroupClass = case_when(isPN & CNAclassification_step1 == "Malignant" ~ "Malignant",
                                isPN & CNAclassification_step1 == "Non-malignant" ~ "PM-NM",
                                isPN & CNAclassification_step1 == "Unresolved" ~ "PM-U",
                                isReference & CNAclassification_step1 == "Non-malignant" ~ "REF-NM",
                                isReference & CNAclassification_step1 == "Unresolved" ~ "REF-U",
                                isONM & CNAclassification_step1 == "Non-malignant" ~ "ONM-NM",
                                isONM & CNAclassification_step1 == "Unresolved" ~ "ONM-U",
                                TRUE ~ "NA")) %>%
  mutate(GroupClass = factor(GroupClass, c("Malignant", "PM-NM", "REF-NM", "ONM-NM", "PM-U", "REF-U", "ONM-U")))

table(cell_classification_data$CNAclassification_step1)
table(cell_classification_data$CellType, cell_classification_data$CNAclassification_step1)
table(cell_classification_data$CellType, cell_classification_data$GroupClass)

table(cell_classification_data$CellType, cell_classification_data$GroupClass) %>% as.matrix() %>%
  write.csv(paste0(TABLES_ROOT, "cell_classification_data_groupclass_vs_prevclass.csv"), row.names = T, col.names = T)

cell_classification_data$Sample <- get_sname(cell_classification_data$CellID, split = "-")

setdiff(meta_data %>% filter(CellType == "Malignant", Sample %in% cell_classification_data$Sample) %>% pull(CellID),
        cell_classification_data %>% filter(CNAclassification_step1 == "Malignant") %>% pull(CellID)) %>%
  length()

meta_data %>% filter(CellType == "Malignant", Sample %in% cell_classification_data$Sample) %>% pull(CellID) %>% length()

setdiff(cell_classification_data %>% filter(CNAclassification_step1 == "Malignant") %>% pull(CellID),
        meta_data %>% filter(CellType == "Malignant", Sample %in% cell_classification_data$Sample) %>% pull(CellID)) %>%
  length()

cell_classification_data %>% filter(CNAclassification_step1 == "Malignant") %>% pull(CellID) %>% length()

setdiff(cell_classification_data %>% filter(CNAclassification_step1 == "Non-malignant") %>% pull(CellID),
        meta_data %>% filter(CNAclassification == "Non-malignant", Sample %in% cell_classification_data$Sample) %>% pull(CellID)) %>%
  length()

#############################################################################################################################################################
# CNA classification STEP 2 - Support classification with CNA correlation
#############################################################################################################################################################

control_panel <- lapply(cna_matrix_list_refcorrected, colnames)
lengths(control_panel)
length(control_panel)
control_panel <- unname(unlist(control_panel))
length(control_panel)
control_panel <- table(control_panel) %>% sort(decreasing = T)
control_panel <- control_panel[control_panel == max(control_panel)] %>% names()
length(control_panel)

i <- 0
cna_cor_per_chr <- lapply(cna_matrix_list_refcorrected, function(cna_matrix) {
  
  i <<- i + 1
  
  print(paste0("Iteration ", i))
  
  if(is.null(cna_matrix))
    return(NULL)
  
  d_sig <- cna_sig_per_chr %>%
    filter(CellID %in% rownames(cna_matrix))
  dim(d_sig)
  
  d_sig <- d_sig %>%
    filter(GroupClass == "Malignant")
  dim(d_sig)
  
  if(nrow(d_sig) == 0) {
    print("No malignant cells were found")
    return(NULL)
  }
  
  cna_events <- d_sig$CHR[d_sig$isRealEvent] %>% unique() %>% as.character()
  tbl <- table(d_sig$CHR[d_sig$isRealEvent]) / sum(table(d_sig$CHR[d_sig$isRealEvent]))
  cna_events <- names(tbl[tbl > .05])
  
  print(paste0(length(cna_events), " high frequency CNA events detected"))
  
  stopifnot(length(cna_events) > 0)
  
  cna_genes <- splitGenes(colnames(cna_matrix))
  cna_genes <- cna_genes[-c(length(cna_genes))]
  names(cna_genes) <- paste0("CHR", names(cna_genes))
  
  cna_genes <- cna_genes[cna_events]
  length(cna_genes)
  
  cna_genes <- unname(unlist(cna_genes, recursive = F))
  
  cp <- setdiff(control_panel, cna_genes)
  
  print(paste0("Adding ", length(cp), " control genes for CNA correlation"))
  
  cna_genes <- c(cna_genes, cp)
  
  d <- meta_data %>%
    filter(CellID %in% d_sig$CellID)
  
  cna0 <- colMeans(cna_matrix[d$CellID, colnames(cna_matrix) %in% cna_genes])
  
  cna_cor <- cor(t(cna_matrix[, colnames(cna_matrix) %in% cna_genes]), cna0)
  
  res <- tibble(CellID = rownames(cna_matrix), CNAcor = cna_cor[, 1])
  
  return(res)
})
cna_cor_per_chr <- do.call(rbind, cna_cor_per_chr)

cell_classification_data <- cell_classification_data %>%
  left_join(cna_cor_per_chr, by = "CellID")

ggplot(cell_classification_data, aes(x = CNAcor, y = GroupClass, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE) +
  scale_fill_viridis_d(name = "Quartiles") +
  xlab("CNA correlation") +
  ylab("") +
  scale_x_continuous(limits = c(-1, 1)) +
  theme_gbm_pvsr(panel.border = element_blank())

#############################################################################################################################################################
# CNA classification STEP 2 - Estimate algorithm performance using mutation calls and set CNAcor threshold
#############################################################################################################################################################

mut_data_path <- "/mut_call/"

snv_per_cell_data <- lapply(list.dirs(mut_data_path), function(dir) {
  sample <- strsplit(dir, "/")
  sample <- sample[[1]][length(sample[[1]])]
  
  print(sample)
  
  f <- paste0(dir, "/", sample, "_SNP.mtx")
  if(!file.exists(f))
    return(NULL)
  
  # Read in the sparse genotype matrix
  snv_matrix <- readMM(f)
  dim(snv_matrix)
  snv_matrix[1:5, 1:5]
  
  # convert the matrix to a dataframe
  snv_matrix <- as.data.frame(as.matrix(snv_matrix))
  
  f <- paste0(dir, "/", sample, "_barcodes.tsv")
  if(!file.exists(f))
    return(NULL)
  
  #read in the cell barcodes output by Cell Ranger
  barcodes <- read.table(f, header = F)
  barcodes$V1 <- paste0(sample, "-", barcodes$V1)
  dim(barcodes)
  
  # Construct the final table to add to the Seurat object
  f <- paste0(dir, "/", sample, "_SNV_loci.txt")
  if(!file.exists(f))
    return(NULL)
  
  snps <- read.table(f, header = F)
  dim(snps)
  snps$V1 <- paste0("CHR", snps$V1)
  
  colnames(snv_matrix) <- barcodes$V1
  
  
  snv_calls <- lapply(1:ncol(snv_matrix), function(i) which(snv_matrix[, i] > 1))
  snv_calls <- unlist(snv_calls, recursive = F)
  snv_calls <- snps[snv_calls, ]
  
  snv_per_cell <- colSums(snv_matrix > 1)
  table(snv_per_cell) %>% sort(decreasing = T)
  
  res <- tibble(CellID = colnames(snv_matrix), Sample = sample, SNVperCell = snv_per_cell)
  
  return(res)
})
snv_per_cell_data <- do.call(rbind, snv_per_cell_data)
dim(snv_per_cell_data)

table(snv_per_cell_data$Sample) %>% sort(decreasing = T)
table(snv_per_cell_data$SNVperCell) %>% sort(decreasing = T)

snv_per_cell_data <- snv_per_cell_data %>%
  filter(CellID %in% cell_classification_data$CellID)
dim(snv_per_cell_data)

d <- cell_classification_data %>%
  left_join(snv_per_cell_data %>%
              select(CellID, SNVperCell),
            by = "CellID") %>%
  filter(!is.na(SNVperCell))
dim(d)

d <- d %>%
  mutate(isSNV = SNVperCell > 0, isSNV_fct = ifelse(isSNV, "SNV detected", "SNV not detected")) %>%
  mutate(isSNV_fct = factor(isSNV_fct, c("SNV detected", "SNV not detected")))

d %>%
  group_by(Sample, GroupClass) %>%
  summarise(n = sum(SNVperCell > 0), N = n(), SNVperCell = n / N) %>%
  ggboxplot(x = "GroupClass", y = "SNVperCell", add = "jitter", xlab = "", ylab = "% cells with detected SNV per sample") +
  scale_y_continuous(labels = percent) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())

d %>%
  group_by(Sample, CNAclassification_step1) %>%
  summarise(n = sum(SNVperCell > 0), N = n(), SNVperCell = n / N) %>%
  ungroup() %>%
  filter(CNAclassification_step1 == "Malignant") %>%
  summarise(Q75 = quantile(SNVperCell, .75, na.rm = T))

d %>%
  group_by(Sample, CNAclassification_step1) %>%
  summarise(n = sum(SNVperCell > 0), N = n(), SNVperCell = n / N) %>%
  filter(CNAclassification_step1 == "Non-malignant") %>%
  arrange(desc(SNVperCell))

d_stats <- d %>%
  filter(!is.na(SNVperCell)) %>%
  group_by(Sample, GroupClass) %>%
  summarise(n = sum(SNVperCell > 0), N = n(), Freq = n / N)

sample_ord <- d_stats %>%
  filter(GroupClass == "Malignant") %>%
  arrange(Freq) %>%
  pull(Sample)

d_stats$Sample <- factor(d_stats$Sample, sample_ord)

d_stats %>%
  filter(GroupClass %in% c("Malignant", "REF-NM", "PM-NM")) %>%
  ggline(x = "Sample", y = "Freq", color = "GroupClass", xlab = "Sample", ylab = "% cells with detected SNV") +
  geom_hline(yintercept = .025, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = .05, linetype = "dashed", color = "black", size = 1) +
  scale_color_brewer(name = "", palette = "Set1") +
  scale_y_continuous(labels = percent) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())

cor(d_stats$Freq[d_stats$GroupClass == "Malignant"], d_stats$Freq[d_stats$GroupClass == "REF-NM"])
cor(d_stats$Freq[d_stats$GroupClass == "Malignant"], d_stats$Freq[d_stats$GroupClass == "PM-NM"])
cor(d_stats$Freq[d_stats$GroupClass == "PM-NM"], d_stats$Freq[d_stats$GroupClass == "REF-NM"])

d_stats <- d %>%
  filter(!is.na(SNVperCell)) %>%
  group_by(Sample, GroupClass) %>%
  summarise(n = sum(SNVperCell > 0), N = n(), Freq = n / N) %>%
  ungroup() %>%
  select(Sample, GroupClass, Freq) %>%
  left_join(d %>%
              filter(!is.na(SNVperCell), isSNV) %>%
              group_by(Sample, GroupClass) %>%
              summarise(n = sum(SNVperCell > 0)) %>%
              group_by(GroupClass) %>%
              mutate(N = sum(n), Weight = n / N) %>%
              ungroup() %>%
              select(Sample, GroupClass, Weight),
            by = c("Sample", "GroupClass"))
d_stats$Weight[is.na(d_stats$Weight)] <- 0
d_stats$Freq[is.na(d_stats$Freq)] <- 0

d_stats %>%
  filter(GroupClass %in% c("PM-NM", "REF-NM"), Weight > .025) %>%
  group_by(GroupClass) %>%
  summarise(Sum = sum(Weight), n = n())

d_stats %>%
  filter(GroupClass %in% c("REF-NM", "PM-NM")) %>%
  mutate(GroupClass = relevel(GroupClass, "REF-NM")) %>%
  ggscatter(x = "Freq", y = "Weight", xlab = "% cells with detected SNV", ylab = "SNV weight [%]") +
  facet_wrap(~GroupClass) +
  geom_vline(xintercept = .05, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = .025, linetype = "dashed", color = "black", size = 1) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent) +
  theme_gbm_pvsr(panel.spacing = unit(1, "lines")) +
  theme(panel.grid.major = element_line())

d_stats %>%
  ungroup() %>%
  filter(GroupClass == "REF-NM") %>%
  summarise(n_025 = sum(Freq < .025), n_05 = sum(Freq < .05), N = n(), Freq_025 = n_025 / N, Freq_05 = n_05 / N)

d_stats %>%
  ungroup() %>%
  filter(GroupClass == "REF-NM") %>%
  summarise(n_025 = sum(Weight < .025), n_05 = sum(Weight < .05), N = n(), Weight_025 = n_025 / N, Weight_05 = n_05 / N)

suspicious_mut_samples <- d_stats %>%
  filter(GroupClass == "REF-NM", Weight > .025) %>%
  pull(Sample)

d <- d %>%
  mutate(isSuspicoius = Sample %in% suspicious_mut_samples,
         isSuspicoius_fct = ifelse(isSuspicoius, "Suspicious sample", "Non-suspicious sample")) %>%
  mutate(isSuspicoius_fct = as.factor(isSuspicoius_fct))

d %>%
  mutate(isSuspicious = Sample %in% suspicious_mut_samples) %>%
  group_by(Sample, GroupClass, isSuspicious) %>%
  summarise(n = sum(SNVperCell> 0), N = n(), SNVperCell = n / N) %>%
  group_by(GroupClass, isSuspicious) %>%
  summarise(Mean = mean(SNVperCell, na.rm = T))

d %>%
  ungroup() %>%
  group_by(Sample, GroupClass, isSNV_fct, isSuspicoius_fct) %>%
  summarise(CNAcor = quantile(CNAcor, .5)) %>%
  ggboxplot(x = "GroupClass", y = "CNAcor", color = "isSNV_fct", add = "jitter",
            xlab = "", ylab = "Median CNA correlation per sample") +
  scale_color_brewer(name = "", palette = "Set1") +
  scale_y_continuous(limits = c(-1, 1)) +
  geom_hline(yintercept = .5, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = .35, linetype = "dashed", color = "black", size = 1) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())

th_data <- d %>%
  ungroup() %>%
  group_by(Sample, GroupClass, isSNV, isSuspicoius) %>%
  summarise(CNAcor = quantile(CNAcor, .5)) %>%
  group_by(GroupClass, isSNV, isSuspicoius) %>%
  summarise(Median = median(CNAcor), Q25 = quantile(CNAcor, .25), Q75 = quantile(CNAcor, .75))

lt <- th_data$Q25[th_data$GroupClass == "PM-NM" & th_data$isSNV == T & th_data$isSuspicoius == F]
ht <- th_data$Q75[th_data$GroupClass == "PM-NM" & th_data$isSNV == T & th_data$isSuspicoius == F]
lt; ht

nf <- MASS::fitdistr(d %>%
                       filter(GroupClass == "REF-NM", isSNV_fct != "SNV detected") %>%
                       pull(CNAcor),
                     "normal")

d_tmp <- rbind(d %>%
                 select(CNAcor, GroupClass, isSNV_fct) %>%
                 mutate(isSNV_fct = as.character(isSNV_fct)),
               tibble(CNAcor = rnorm(n = 100000, mean = nf$estimate[1], sd = nf$estimate[2]),
                      GroupClass = "REF-NM", isSNV_fct = "Classifier"),
               tibble(CNAcor = rnorm(n = 100000, mean = nf$estimate[1], sd = nf$estimate[2]),
                      GroupClass = "PM-NM", isSNV_fct = "Classifier")) %>%
  mutate(isSNV_fct = factor(isSNV_fct, c("SNV not detected", "SNV detected", "Classifier")))

d_tmp %>%
  filter(GroupClass %in% c("PM-NM", "REF-NM")) %>%
  mutate(GroupClass = relevel(GroupClass, "REF-NM")) %>%
  ggplot(aes(x = CNAcor, y = after_stat(ncount), color = isSNV_fct, linetype = isSNV_fct)) +
  facet_wrap(~GroupClass) +
  geom_freqpoly(binwidth = .01, show.legend = c(color = T, linetype = F), size = 1) +
  scale_color_manual(name = "", values = c("SNV not detected" = "dodgerblue", "SNV detected" = "red", "Classifier" = "black")) +
  scale_linetype_manual(values = c("SNV not detected" = "solid", "SNV detected" = "solid", "Classifier" = "dashed")) +
  geom_vline(xintercept = nf$estimate[1] + 1.96*nf$estimate[2], linetype = "dashed", color = "black", size = 1) +
  scale_x_continuous(limits = c(-1, 1)) +
  xlab("CNA correlation") +
  ylab("Count (scaled to 1)") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_gbm_pvsr() +
  theme(panel.grid.major.x = element_line())

lt <- nf$estimate[1] + 1.96*nf$estimate[2]
ht <- .5

# Estimate lower threshold for discarding cells
cor_th_est <- lapply(seq(0, 1, .05), function(cor_th) {
  
  d_tmp <- d %>%
    filter(!isSuspicoius)
  
  d_tmp <- d_tmp %>%
    mutate(CNAclassification_step2 = case_when(CNAclassification_step1 == "Malignant" & CNAcor < ht ~ "Unresolved",
                                               CNAclassification_step1 == "Non-malignant" & CNAcor > cor_th ~ "Unresolved",
                                               TRUE ~ CNAclassification_step1))
  table(d_tmp$CNAclassification_step1, d_tmp$CNAclassification_step2)
  table(d_tmp$CNAclassification_step2, d_tmp$GroupClass)
  
  d_tmp$GroupClass2 <- as.character(d_tmp$GroupClass)
  
  d_tmp <- d_tmp %>%
    mutate(GroupClass2 = case_when(CNAclassification_step2 == "Unresolved" & GroupClass2 == "Malignant" ~ "PM-U",
                                   CNAclassification_step2 == "Unresolved" & GroupClass2 == "PM-NM" ~ "PM-U",
                                   CNAclassification_step2 == "Unresolved" & GroupClass2 == "REF-NM" ~ "REF-U",
                                   CNAclassification_step2 == "Unresolved" & GroupClass2 == "ONM-NM" ~ "ONM-U",
                                   TRUE ~ GroupClass2)) %>%
    mutate(GroupClass2 = factor(GroupClass2, levels(GroupClass)))
  table(d_tmp$GroupClass, d_tmp$GroupClass2)
  
  d_tmp_stats1 <- d_tmp %>%
    group_by(GroupClass) %>%
    summarise(n = sum(CNAclassification_step2 == "Unresolved"), N = n(), dFiscarded = n / N)
  d_tmp_stats1
  
  d_tmp_stats2 <- d_tmp %>%
    group_by(GroupClass2) %>%
    summarise(n = sum(SNVperCell > 0), N = n(), fMissclassified = n / N, CNAcor = median(CNAcor))
  d_tmp_stats2
  
  d_tmp_stats3 <- d_tmp %>%
    filter(!isSNV) %>%
    group_by(GroupClass2) %>%
    summarise(CNAcor = median(CNAcor))
  d_tmp_stats3
  
  d_tmp_stats4 <- d_tmp %>%
    filter(isSNV) %>%
    group_by(GroupClass2) %>%
    summarise(CNAcor = median(CNAcor))
  d_tmp_stats4
  
  res <- tibble(CorTH = cor_th, GroupClass2 = d_tmp_stats1$GroupClass, dFiscarded = d_tmp_stats1$dFiscarded, fMissclassified = d_tmp_stats2$fMissclassified, CNAcor_woSNV = d_tmp_stats3$CNAcor, CNAcor_wSNV = d_tmp_stats4$CNAcor)
  
  return(res)  
})
cor_th_est <- do.call(rbind, cor_th_est)

p1 <- cor_th_est %>%
  filter(GroupClass2 %in% c("PM-NM", "REF-NM")) %>%
  ggline(x = "CorTH", y = "fMissclassified", color = "GroupClass2", xlab = "CNA correlation threshold", ylab = "Percent cells missclassified") +
  geom_vline(xintercept = .3, linetype = "dashed", color = "black") +
  scale_color_brewer(name = "", palette = "Set1") +
  scale_y_continuous(labels = percent) +
  scale_x_discrete(breaks = seq(0, 1, .1)) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())
p2 <- cor_th_est %>%
  filter(GroupClass2 %in% c("PM-NM", "REF-NM")) %>%
  ggline(x = "CorTH", y = "dFiscarded", color = "GroupClass2", xlab = "CNA correlation threshold", ylab = "Percent cells discarded") +
  geom_vline(xintercept = .3, linetype = "dashed", color = "black") +
  scale_color_brewer(name = "", palette = "Set1") +
  scale_x_discrete(breaks = seq(0, 1, .1)) +
  scale_y_continuous(labels = percent) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())
p3 <- cor_th_est %>%
  filter(GroupClass2 %in% c("PM-NM", "REF-NM")) %>%
  ggline(x = "CorTH", y = "CNAcor_woSNV", color = "GroupClass2", xlab = "CNA correlation threshold",
         ylab = "Median CNA correlation\n(no SNV)") +
  scale_color_brewer(name = "", palette = "Set1") +
  scale_x_discrete(breaks = seq(0, 1, .1)) +
  scale_y_continuous(limits = c(-.2, .5)) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())
p4 <- cor_th_est %>%
  filter(GroupClass2 %in% c("PM-NM", "REF-NM")) %>%
  ggline(x = "CorTH", y = "CNAcor_wSNV", color = "GroupClass2", xlab = "CNA correlation threshold",
         ylab = "Median CNA correlation (with SNV)") +
  scale_color_brewer(name = "", palette = "Set1") +
  scale_x_discrete(breaks = seq(0, 1, .1)) +
  scale_y_continuous(limits = c(-.2, .5)) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

p1 + p2 + p3 + p4

p1
p2

cell_classification_data <- cell_classification_data %>%
  left_join(snv_per_cell_data %>%
              select(CellID, SNVperCell),
            by = "CellID")

cell_classification_data <- cell_classification_data %>%
  mutate(CNAclassification_step2 = case_when(CNAclassification_step1 == "Malignant" & CNAcor < ht ~ "Unresolved",
                                             CNAclassification_step1 == "Non-malignant" & CNAcor > lt ~ "Unresolved",
                                             TRUE ~ CNAclassification_step1))
table(cell_classification_data$CNAclassification_step1, cell_classification_data$CNAclassification_step2)
table(cell_classification_data$CNAclassification_step2, cell_classification_data$GroupClass)

tbl <- table(cell_classification_data$CNAclassification_step1, cell_classification_data$CNAclassification_step2)
tbl[2, 3] / (tbl[2, 2] + tbl[2, 3])

cell_classification_data$GroupClass2 <- as.character(cell_classification_data$GroupClass)

cell_classification_data <- cell_classification_data %>%
  mutate(GroupClass2 = case_when(CNAclassification_step2 == "Unresolved" & GroupClass2 == "Malignant" ~ "PM-U",
                                 CNAclassification_step2 == "Unresolved" & GroupClass2 == "PM-NM" ~ "PM-U",
                                 CNAclassification_step2 == "Unresolved" & GroupClass2 == "REF-NM" ~ "REF-U",
                                 CNAclassification_step2 == "Unresolved" & GroupClass2 == "ONM-NM" ~ "ONM-U",
                                 TRUE ~ GroupClass2)) %>%
  mutate(GroupClass2 = factor(GroupClass2, levels(GroupClass)))
table(cell_classification_data$GroupClass, cell_classification_data$GroupClass2)

cell_classification_data %>%
  ungroup() %>%
  filter(!is.na(SNVperCell)) %>%
  mutate(isSNV = ifelse(SNVperCell > 0, "SNV detected", "SNV not detected")) %>%
  mutate(isSNV = factor(isSNV, c("SNV detected", "SNV not detected"))) %>%
  mutate(isSuspicoius = Sample %in% suspicious_mut_samples,
         isSuspicoius_fct = ifelse(isSuspicoius, "Suspicious sample", "Non-suspicious sample")) %>%
  mutate(isSuspicoius_fct = as.factor(isSuspicoius_fct)) %>%
  group_by(Sample, GroupClass2, isSNV, isSuspicoius_fct) %>%
  summarise(CNAcor = quantile(CNAcor, .5)) %>%
  ggboxplot(x = "GroupClass2", y = "CNAcor", color = "isSNV", add = "jitter",
            xlab = "", ylab = "Median CNA correlation per sample") +
  scale_color_brewer(name = "", palette = "Set1") +
  scale_y_continuous(limits = c(-1, 1)) +
  geom_hline(yintercept = .5, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = .35, linetype = "dashed", color = "black", size = 1) +
  theme_gbm_pvsr(panel.border = element_blank(), axis.text.x = element_text(size = 16, angle = 90)) +
  theme(panel.grid.major = element_line())

cell_classification_data %>% 
  mutate(isSNV = SNVperCell > 0) %>% 
  select(GroupClass2, isSNV) %>%
  table()

#############################################################################################################################################################
# CNA classification STEP 3 - Classification confidence levels
#############################################################################################################################################################

cna_sig_per_chr <- cna_sig_per_chr %>%
  left_join(cell_classification_data %>%
              select(CellID, CNAclassification_step2, GroupClass2),
            by = "CellID")

cna_sig_per_chr_stats <- cna_sig_per_chr %>%
  group_by(Sample, CHR, CNAclassification_step2, GroupClass2) %>%
  summarise(n = sum(isRealEvent), N = n(), Freq = n / N)

cna_sig_per_chr_stats %>%
  filter(CNAclassification_step2 == "Malignant") %>%
  ggboxplot(x = "CHR", y = "Freq", add = "jitter")

cell_classification_data <- cell_classification_data %>% 
  mutate(isSNV = SNVperCell > 0)

n_sample_events <- cna_sig_per_chr %>%
  group_by(Sample, CHR) %>%
  summarise(isEvent = sum(isRealEvent) > 0) %>%
  group_by(Sample) %>%
  summarise(nSampleEvents = sum(isEvent))

cna_sig_per_chr %>%
  group_by(Sample, CHR, GroupClass2) %>%
  summarise(n = sum(isRealEvent), N = n(), Freq = n / N) %>%
  filter(Freq > 0) %>%
  View()

cell_classification_data <- cell_classification_data %>%
  left_join(cna_sig_per_chr %>%
              group_by(CellID) %>%
              summarise(nEvents = sum(isRealEvent)),
            by = "CellID")

cell_classification_data  <- cell_classification_data %>%
  left_join(cna_sig_per_chr %>%
              group_by(CellID) %>%
              summarise(nSuspicious = sum(isSuspiciousEvent)),
            by = "CellID")

cell_classification_data <- cell_classification_data %>%
  left_join(n_sample_events, by = "Sample")

cell_classification_data$eventFreq <- cell_classification_data$nEvents / cell_classification_data$nSampleEvents

cell_classification_data <- cell_classification_data %>%
  left_join(cna_sig_per_chr %>%
              group_by(CellID) %>%
              summarise(isCHR7 = isRealEvent[CHR == "CHR7"], isCHR10 = isRealEvent[CHR == "CHR10"]),
            by = "CellID")

cell_classification_data %>%
  ggboxplot(x = "GroupClass2", y = "eventFreq", add = "jitter")

cell_classification_data %>%
  filter(GroupClass2 == "PM-U") %>%
  ggplot(aes(x = CNAcor, y = as.factor(nSuspicious), fill = factor(stat(quantile)))) +
  facet_wrap(~isSNV) +
  stat_density_ridges(geom = "density_ridges_gradient", color = "black", show.legend = c(fill = F),
                      calc_ecdf = TRUE, quantiles = 4, quantile_lines = TRUE) +
  geom_vline(xintercept = .5, linetype = "dashed", color = "black", size = 1) +
  scale_fill_viridis_d(name = "Quartiles") +
  xlab("CNA correlation") +
  ylab("# suspicious events") +
  theme_gbm_pvsr()

# Return from the dead
cell_classification_data <- cell_classification_data %>%
  mutate(CNAclassification_step3 = case_when(GroupClass2 == "PM-U" & nEvents > 0 ~ "Malignant",
                                             GroupClass2 == "PM-U" & nSuspicious > 0 & CNAcor > .5 ~ "Malignant",
                                             GroupClass2 == "PM-U" & nSuspicious > 0 & isSNV ~ "Malignant",
                                             GroupClass2 == "PM-U" & nSuspicious == 0 & CNAcor < .5 ~ "Non-malignant",
                                             CNAclassification_step2 == "Non-malignant" & isSNV ~ "Unresolved",
                                             TRUE ~ CNAclassification_step2))
table(cell_classification_data$CNAclassification_step3)

cell_classification_data <- cell_classification_data %>%
  mutate(GroupClass3 = as.character(GroupClass2)) %>%
  mutate(GroupClass3 = case_when(CNAclassification_step3 == "Malignant" & GroupClass3 == "PM-U" ~ "Malignant",
                                 CNAclassification_step3 == "Non-malignant" & GroupClass3 == "PM-U" ~ "PM-NM",
                                 CNAclassification_step3 == "Unresolved" & GroupClass3 == "PM-NM" ~ "PM-U",
                                 CNAclassification_step3 == "Unresolved" & GroupClass3 == "REF-NM" ~ "REF-U",
                                 CNAclassification_step3 == "Unresolved" & GroupClass3 == "ONM-NM" ~ "ONM-U",
                                 TRUE ~ GroupClass3)) %>%
  mutate(GroupClass3 = factor(GroupClass3, levels(GroupClass2)))
table(cell_classification_data$CNAclassification_step3, cell_classification_data$GroupClass3)

cell_classification_data <- cell_classification_data %>%
  mutate(ClassConf = case_when(CNAclassification_step3 == "Malignant" & isSNV ~ "1",
                               CNAclassification_step3 == "Malignant" & eventFreq >= .5 & (isCHR7 | isCHR10) ~ "1",
                               CNAclassification_step3 == "Malignant" & eventFreq >= .5 & (!isCHR7 & !isCHR10) ~ "2",
                               CNAclassification_step3 == "Malignant" & eventFreq < .5 & (isCHR7 | isCHR10) ~ "2",
                               CNAclassification_step3 == "Malignant" & nEvents > 0 ~ "3",
                               CNAclassification_step3 == "Malignant" & nSuspicious > 0 & CNAcor > .5 ~ "3",
                               CNAclassification_step3 == "Malignant" & nSuspicious > 0 & isSNV ~ "3",
                               CNAclassification_step3 == "Non-malignant" & CNAcor < .1 ~ "1",
                               CNAclassification_step3 == "Non-malignant" & CNAcor < .25 ~ "2",
                               CNAclassification_step3 == "Non-malignant" & CNAcor < .5 ~ "3",
                               CNAclassification_step3 == "Unresolved" ~ "4",
                               # To identify leaky conditions
                               TRUE ~ "5"))
table(cell_classification_data$CNAclassification_step3, cell_classification_data$ClassConf)
table(cell_classification_data$GroupClass3, cell_classification_data$ClassConf)

# Penalize samples with no WES CNA calls by forcing max conf class to 2
cc_vec <- sapply(1:nrow(cell_classification_data), function(i) ifelse(cell_classification_data$Sample[i] %in% no_call_samples, max(cell_classification_data$ClassConf[i], 2), cell_classification_data$ClassConf[i]))
table(cc_vec)
table(cell_classification_data$ClassConf)

cell_classification_data$ClassConf <- cc_vec

colSums(table(cell_classification_data$CNAclassification_step3, cell_classification_data$ClassConf))

table(cell_classification_data$GroupClass3, cell_classification_data$nEvents)
table(cell_classification_data$GroupClass3, cell_classification_data$nSuspicious)

table(cell_classification_data$CellType, cell_classification_data$ClassConf)
table(cell_classification_data$CellType, cell_classification_data$GroupClass3)

table(cell_classification_data$GroupClass3, cell_classification_data$ClassConf) %>%
  write.csv(paste0(TABLES_ROOT, "class_conf.csv"), row.names = T, col.names = T)

d <- cell_classification_data %>%
  left_join(meta_data %>%
              select(CellID, Patient_factor, Timepoint, ID_factor),
            by = "CellID") %>%
  filter(Timepoint != "T3") %>%
  group_by(ID_factor, Patient_factor, Timepoint, CNAclassification_step3, ClassConf) %>%
  summarise(n = n())

d %>%
  filter(CNAclassification_step3 == "Malignant") %>%
  ggplot(aes(x = Timepoint, y = n, fill = ClassConf)) +
  facet_wrap(~Patient_factor, ncol = 12) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", size = 1) +
  scale_fill_brewer("Confidence level", palette = "Set1") +
  scale_y_log10() +
  xlab("") +
  ylab("# cells [log10]") +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90), panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

d %>%
  filter(CNAclassification_step3 == "Non-malignant") %>%
  ggplot(aes(x = Timepoint, y = n, fill = ClassConf)) +
  facet_wrap(~Patient_factor, ncol = 12) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", size = 1) +
  scale_fill_brewer("Confidence level", palette = "Set1") +
  scale_y_log10() +
  xlab("") +
  ylab("# cells [log10]") +
  theme_gbm_pvsr(axis.text.x = element_text(size = 16, angle = 90), panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

tmp <- acast(d %>% filter(CNAclassification_step3 == "Malignant"), formula = ID_factor ~ ClassConf, value.var = "n")
tmp[is.na(tmp)] <- 0

colMeans(tmp)
apply(tmp, 2, median)
apply(tmp, 2, quantile)

tmp2 <- tibble(ID = rownames(tmp), n = tmp[, 1] + tmp[, 2])

tmp2$Patient <- sapply(strsplit(tmp2$ID, split = "T"), function(x) x[[1]][1])
tmp2$Timepoint <- sapply(strsplit(tmp2$ID, split = "T"), function(x) x[[2]][1])
tmp2$Timepoint <- paste0("T", tmp2$Timepoint)

tmp2 <- tmp2 %>%
  group_by(Patient) %>%
  filter("T1" %in% Timepoint & "T2" %in% Timepoint)
dim(tmp2)

tmp2_stats <- tmp2 %>%
  group_by(Patient) %>%
  summarise(has100 = n[Timepoint == "T1"] >= 100 & n[Timepoint == "T2"] >= 100,
            has50 = n[Timepoint == "T1"] >= 50 & n[Timepoint == "T2"] >= 50)

sum(tmp2_stats$has100) / nrow(tmp2_stats)
sum(tmp2_stats$has50) / nrow(tmp2_stats)

cna_sig_per_chr <- cna_sig_per_chr %>%
  left_join(cell_classification_data %>%
              select(CellID, CNAclassification_step3, GroupClass3, ClassConf),
            by = "CellID")

saveRDS(cell_classification_data, paste0(DATA_ROOT, "cell_classification_data.RDS"))
saveRDS(cna_sig_per_chr, paste0(DATA_ROOT, "cna_sig_per_chr.RDS"))

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Classify the non-malignant cells
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

nmdata <- meta_data %>%
  filter(CellID %in% cell_classification_data$CellID) %>%
  select(CellID, orig.ident, Sample, Patient, Patient_factor, Timepoint, Timepoint_factor, ID, ID_extended, ID_factor, ID_extended_factor, Origin, Lab, Pathology, IDHstatus, MGMTstatus, SurgicalIntervalmonth, OSmonth, PatientStatus, CellType_old)

nmdata <- nmdata %>%
  left_join(cell_classification_data %>%
              select(CellID, nEvents, nSuspicious, CNAcor, CNAclassification, GroupClass, ClassConf),
            by = "CellID")

nmdata <- nmdata %>%
  filter(CNAclassification == "Non-malignant", ClassConf %in% c("1", "2", "3"))
dim(nmdata)
table(nmdata$GroupClass)
table(nmdata$ClassConf)

nmdata <- as.data.frame(nmdata)
rownames(nmdata) <- nmdata$CellID

table(nmdata$CNAclassification_step3)
which(duplicated(nmdata$CellID))
table(nmdata$CellType_old)

#############################################################################################################################################################
# Select genes
#############################################################################################################################################################

selected_genes <- select_genes(umi_data_list = umi_data_all, md = nmdata, exp_th = 1, freq_th2 = 5, verbose = T, plot = T)
length(selected_genes)

junk_genes <- selected_genes[grep("\\.", selected_genes)]
length(junk_genes)

selected_genes <- setdiff(selected_genes, junk_genes)
length(selected_genes)

gc()

#############################################################################################################################################################
# Generate UMI matrix for seurat dataset
#############################################################################################################################################################

table(nmdata$Sample) %>% sort(decreasing = F)

umi_data <- lapply(unique(nmdata$Sample), function(sname) {
  m <- umi_data_all[[sname]]
  m <- m[selected_genes, nmdata$CellID[nmdata$Sample == sname]]
  m
})
umi_data <- do.call(cbind, umi_data)
gc()
dim(umi_data)

umi_data <- umi_data[, nmdata$CellID]
dim(umi_data)

saveRDS(object = umi_data, file = paste0(DATA_ROOT, "umi_data_nm.RDS"))

rm(umi_data_all)
gc()

#############################################################################################################################################################
# Create a top-level Seurat object
#############################################################################################################################################################

nmsds <- new_seurat_object(m = umi_data, md = nmdata, project = "GBM-CARE_integrated_analysis_nonmalignant",
                           nfeatures = 3000, metric = "correlation", verbose = T, npcs = 100)
dim(nmsds)

DimPlot(nmsds, reduction = "umap", group.by = "Patient") +
  DimPlot(nmsds, reduction = "umap", group.by = "Timepoint") +
  DimPlot(nmsds, reduction = "umap", group.by = "GroupClass") +
  DimPlot(nmsds, reduction = "umap", group.by = "CellType_old")

saveRDS(object = nmsds, file = paste0(DATA_ROOT, "nmsds.RDS"))

#############################################################################################################################################################
# Cell type annotation
#############################################################################################################################################################

sigs_list <- load_brain_signatures()

sigs_list <- lapply(sigs_list, function(sig) sig[sig %in% rownames(nmsds)])
lengths(sigs_list) %>% sort()
sigs_list <- sigs_list[-c(grep("G1S|G2M", names(sigs_list)))]
sigs_list <- sigs_list[names(sigs_list) %ni% c("IPC.MGE.1-nowak17", "IPC.MGE.2-nowak17")]
length(sigs_list)

obj <- AddModuleScore(object = nmsds, features = sigs_list)

scores <- obj@meta.data[, grep("Cluster", colnames(obj@meta.data))]
colnames(scores) <- names(sigs_list)
dim(scores)

rm(obj)
gc()

saveRDS(scores, paste0(DATA_ROOT, "non_malignant_scores_inter_tumor.RDS"))

nmsds$maxScore <- apply(scores, 1, max)
nmsds$maxSig <- apply(scores, 1, function(x) names(sigs_list)[which.max(x)])

gghistogram(data = nmsds@meta.data, x = "maxScore", color = "black", bins = 50)

null_dist <- generate_null_dist(umi_data_list = umi_data_all, md = nmdata, sigs = sigs_list, genes_subset = rownames(nmsds), verbose = T)

d <- rbind(tibble(Score = nmsds$maxScore, DataType = "Original"),
           tibble(Score = apply(null_dist, 1, max), DataType = "Shuffled"))

ggplot(d, aes(x = Score, y = after_stat(ncount))) +
  facet_wrap(~DataType, nrow = 2) +
  geom_histogram(aes(fill = DataType), color = "black", bins = 50) +
  geom_vline(xintercept = quantile(d$Score[d$DataType == "Shuffled"], .999), color = "black", size = 1, linetype = "dashed") +
  scale_fill_brewer(name = "", palette = "Set1") +
  xlab("Score") +
  ylab("Count (scaled to 1)") +
  theme_gbm_pvsr() +
  theme(panel.grid.major = element_line())

saveRDS(null_dist, paste0(DATA_ROOT, "non_malignant_scores_null_dist.RDS"))

scores_nd <- lapply(colnames(null_dist), function(mp) {
  
  x <- null_dist %>% pull(mp)
  
  fit <- MASS::fitdistr(x, "normal")
  class(fit)
  
  para <- fit$estimate
  
  tibble(MP = mp, Mean = para[1], SD = para[2])  
})
scores_nd <- do.call(rbind, scores_nd)

mean_vec <- setNames(scores_nd$Mean, scores_nd$MP)
sd_vec <- setNames(scores_nd$SD, scores_nd$MP)

ct_data <- melt(data = as_tibble(scores, rownames = "CellID"))
ct_data <- as_tibble(ct_data)
colnames(ct_data) <- c("CellID", "Program", "Score")
ct_data$Program <- as.character(ct_data$Program)
table(ct_data$Program)

ct_data <- ct_data %>%
  mutate(Score_z = (Score - mean_vec[Program]) / sd_vec[Program])

ct_data <- ct_data %>%
  mutate(p.val = pnorm(Score_z, lower.tail = F))

ct_data <- ct_data %>%
  group_by(CellID) %>%
  mutate(p.adj = p.adjust(p.val, "holm"))

ct_data$p.sig <- ct_data$p.adj < .05

ct_data_stats <- ct_data %>%
  group_by(CellID) %>%
  summarise(nsig = sum(p.sig))
table(ct_data_stats$nsig)

ct_data <- ct_data %>%
  ungroup() %>%
  filter(p.sig) %>%
  group_by(CellID) %>%
  arrange(desc(Score)) %>%
  filter(!duplicated(CellID)) %>%
  ungroup()
nrow(ct_data)
table(ct_data$Program) %>% sort()

ct_data$CellType <- as.character(ct_data$Program)
ct_data$CellType[grep("OPC", ct_data$CellType)] <- "OPC"
ct_data$CellType[grep("OC", ct_data$CellType)] <- "Oligodendrocyte"
ct_data$CellType[grep("END", ct_data$CellType)] <- "Endothel"
ct_data$CellType[grep("MIC", ct_data$CellType)] <- "Macrophage"
ct_data$CellType[grep("MAC", ct_data$CellType)] <- "Macrophage"
ct_data$CellType[grep("PVMM", ct_data$CellType)] <- "Macrophage"
ct_data$CellType[grep("PER", ct_data$CellType)] <- "Pericyte"
ct_data$CellType[grep("MUR", ct_data$CellType)] <- "Pericyte"
ct_data$CellType[grep("AC", ct_data$CellType)] <- "Astrocyte"
ct_data$CellType[grep("NEU.EX|NEU.nEX", ct_data$CellType)] <- "Excitatory neuron"
ct_data$CellType[grep("NEU.IN|NEU.nIN", ct_data$CellType)] <- "Inhibitory neuron"
ct_data$CellType[grep("NRGN|GLYC|Cck|GABAergic", ct_data$CellType)] <- "Other neuron"
ct_data$CellType[grep("RG|BG", ct_data$CellType)] <- "RG"
ct_data$CellType[grep("IPC", ct_data$CellType)] <- "IPC"
ct_data$CellType[grep("NEU.|NEU-", ct_data$CellType)] <- "Excitatory neuron"
ct_data$CellType[grep("Choroid|EPN", ct_data$CellType)] <- "Ependyma/Choroid"
ct_data$CellType[grep("ParsTuber", ct_data$CellType)] <- "ParsTuber"
ct_data$CellType[grep("VLMC", ct_data$CellType)] <- "VLMC"

ct_tbl <- table(ct_data$CellType)

sort(ct_tbl, decreasing = T)

ct_data$CellType[ct_data$CellType %in% names(ct_tbl[ct_tbl <= 1200]) & ct_data$CellType != "Bcell"] <- "Other normal"

table(ct_data$CellType) %>% sort(decreasing = T)

unres_cells <- ct_data_stats %>%
  filter(nsig == 0)

celltype_vec <- c(setNames(ct_data$CellType, ct_data$CellID),
                  setNames(rep("Unresolved", nrow(unres_cells)), unres_cells$CellID))
length(celltype_vec)
which(is.na(celltype_vec))

Idents(nmsds) <- factor(celltype_vec[colnames(nmsds)], c("Astrocyte", "Oligodendrocyte", "OPC", "Excitatory neuron", "Inhibitory neuron",
                                                         "Endothel", "Pericyte", "Macrophage", "Tcell", "Bcell", "RG",
                                                         "Other neuron", "Other normal", "Unresolved"))

nmsds$CellType <- Idents(nmsds)

nmdata <- as.data.frame(nmsds@meta.data)

nmdata$UMAP1 <- Embeddings(nmsds, "umap")[, 1]
nmdata$UMAP2 <- Embeddings(nmsds, "umap")[, 2]

DimPlot(nmsds, reduction = "umap", group.by = "CellType") +
  ggtitle("Cell type") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "right")

DimPlot(subset(nmsds, CellType == "OPC"), reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Cell type") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "right")

DimPlot(nmsds, reduction = "umap", group.by = "isCC") +
  ggtitle("Is cycling") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "right")

DimPlot(nmsds, reduction = "umap", group.by = "Lab") +
  ggtitle("Lab") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "right")

DimPlot(nmsds, reduction = "umap", group.by = "IDHstatus") +
  ggtitle("IDH status") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(panel.border = element_blank(), legend.position = "right")

DimPlot(nmsds, reduction = "umap", group.by = "Lab") +
  DimPlot(nmsds, reduction = "umap", group.by = "CellType")

saveRDS(object = nmdata, file = paste0(DATA_ROOT, "nmdata.RDS"))

saveRDS(object = nmsds, file = paste0(DATA_ROOT, "nmsds.RDS"))

####################################################################################################################################
# Clean environment
####################################################################################################################################

rm(cc_null_dist)
rm(cc_null_dist_stats)
rm(cc_scores)
rm(cc_sigs)
rm(ct_tbl)
rm(ct_data)
rm(ct_data_stats)
rm(d)
rm(es_data)
rm(null_dist)
rm(null_dist_stats)
rm(scores)
rm(sigs_list)
rm(unres_cells)
rm(nmsds)
gc()

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Assign cell type and run doublet detection
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

meta_data$CellType_old <- meta_data$CellType

meta_data$CellType <- meta_data$CNAclassification
table(meta_data$CellType)

celltype_vec <- setNames(as.character(nmdata$CellType), nmdata$CellID)

rownames(meta_data) <- meta_data$CellID
meta_data[names(celltype_vec), "CellType"] <- celltype_vec
table(meta_data$CNAclassification, meta_data$CellType)
table(meta_data$IDHstatus, meta_data$CellType)

doublet_det_res <- run_scDblFinder(sample_path = SDS_SAMPLE_PATH, md = meta_data, samples = unique(meta_data$Sample), verbose = T)
table(doublet_det_res$scDblFinder.cluster, doublet_det_res$scDblFinder.class)
table(doublet_det_res$scDblFinder.cluster, doublet_det_res$scDblFinder.class) / rowSums(table(doublet_det_res$scDblFinder.cluster, doublet_det_res$scDblFinder.class)) * 100
colSums(table(doublet_det_res$scDblFinder.cluster, doublet_det_res$scDblFinder.class)) / nrow(doublet_det_res)

saveRDS(doublet_det_res, paste0(DATA_ROOT, "doublet_det_res_scdblfinder_intra_tumor.RDS"))

ggboxplot(data = doublet_det_res, x = "scDblFinder.cluster", y = "scDblFinder.score") +
  theme(axis.text.x = element_text(angle = 90))
ggboxplot(data = doublet_det_res, x = "scDblFinder.class", y = "scDblFinder.score") +
  theme(axis.text.x = element_text(angle = 90))

doublets_vec <- setNames(doublet_det_res$scDblFinder.class, doublet_det_res$CellID)

meta_data$doubletClass <- doublets_vec[meta_data$CellID]

ggscatter(data = meta_data, x = "UMAP1", y = "UMAP2", color = "doubletClass", size = .1)

nmdata$doubletClass <- doublets_vec[nmdata$CellID]

ggscatter(data = nmdata, x = "UMAP1", y = "UMAP2", color = "CellType", size = .1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggscatter(data = nmdata, x = "UMAP1", y = "UMAP2", color = "doubletClass", size = .1) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

meta_data$CellType[meta_data$doubletClass == "doublet"] <- "Doublet"
table(meta_data$CellType)

meta_data %>%
  ungroup() %>%
  mutate(isMalignant = CellType == "Malignant") %>%
  group_by(Sample) %>%
  summarise(n = sum(isMalignant), IDHstatus = first(IDHstatus), Patient = first(Patient), Timepoint = first(Timepoint)) %>%
  arrange(n) %>%
  View()

nmdata$CellType[nmdata$doubletClass == "doublet"] <- "Doublet"
table(nmdata$CellType)

saveRDS(meta_data, paste0(DATA_ROOT, "meta_data.RDS"))
saveRDS(nmdata, paste0(DATA_ROOT, "nmdata.RDS"))

ggscatter(data = nmdata, x = "UMAP1", y = "UMAP2", color = "CellType", size = .1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_gbm_pvsr(panel.border = element_blank())

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Derive NMF meta-programs
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

####################################################################################################################################
# Prep NMF matrices
####################################################################################################################################

md <- readRDS(paste0(DATA_ROOT, "final_cell_classification_data.RDS"))

md <- md %>%
  filter(isValidCell, ClassConf %in% c("1", "2")) %>%
  filter(CellType %ni% c("Unresolved", "LQ", "Other normal")) %>%
  mutate(CellType = as.character(CellType))
dim(md)
table(md$CellType)

md$CellType[md$CellType == "Excitatory neuron"] <- "ExN"
md$CellType[md$CellType == "Inhibitory neuron"] <- "InN"
md$CellType[md$CellType == "Other neuron"] <- "OtN"
table(md$CellType)

celltypes <- unique(md$CellType)
length(celltypes)

samples <- unique(md$Sample)
length(samples)

nmf_matrices_path <- "nmf/nmf_mats/"

start_time <- Sys.time()
for(ct in celltypes) {
  
  save_dir <- paste0(nmf_matrices_path, ct, "/")
  
  if(!dir.exists(paths = save_dir)) {
    print(paste0("Creating directory for ", ct))
    dir.create(path = save_dir)
  }
  
  print("Preparing expression matrices for NMF")
  
  for(i in 1:length(samples)) {
    
    gc()
    
    sname <- samples[i]
    
    print(paste0("******************************* ", sname, " - start, i = ", i, " *******************************"))
    
    d <- md %>%
      filter(Sample == sname, CellType == ct)
    
    if(nrow(d) < 10) {
      print(paste0("Less than 10 ", ct, " cells found for sample ", sname))
      next
    }
    
    m <- umi_data_all[[sname]]
    m <- m[, d$CellID]
    m <- umi2upm(m)
    dim(m)
    
    print(paste0("Found ", ncol(m), " cells"))
    
    rm <- log2(rowMeans(m) + 1)
    
    genes <- names(rm[rm > 4])
    length(genes)
    
    print(paste0("Found ", length(genes), " highly expressed genes"))
    
    genes <- genes[genes %in% rownames(m)]
    m <- m[genes, ]
    m <- log2(m / 10 + 1)
    m <- apply(m, 1, function(x) x - mean(x))
    m <- t(m)
    dim(m)
    
    m[m < 0] <- 0
    dim(m)
    
    m <- m[rowSums(m) > 0, ]
    m <- m[, colSums(m) > 0]
    dim(m)
    
    print("Saving")
    st <- Sys.time()
    saveRDS(object = m, file = paste0(save_dir, sname, ".RDS"))
    et <- Sys.time()
    print(et - st)
    
    print(paste0("******************************* ", sname, " - end, i = ", i, " *******************************"))
  }
}
end_time <- Sys.time()
end_time - start_time

####################################################################################################################################
# Run NMF on each matrix and save output to disk
#
# This code will run NMF for each cell type in each sample. This can run in Rstudio but will take forever, so better to take this
# code snippet and run externally of a high power computation cluster.
#
####################################################################################################################################

nmf_res_path <- "nmf/nmf_res/"

start_time <- Sys.time()
for(ct in celltypes) {
  
  read_dir <- paste0(nmf_matrices_path, ct, "/")
  save_dir <- paste0(nmf_res_path, ct, "/")
  
  if(!dir.exists(paths = save_dir)) {
    print(paste0("Creating directory for ", ct))
    dir.create(path = save_dir)
  }
  
  for(i in 1:length(samples)) {
    
    gc()
    
    sname <- samples[i]
    
    print(paste0("******************************* ", sname, " - start, i = ", i, " *******************************"))
    
    m <- readRDS(file = paste0(read_dir, sname, ".RDS"))
    
    nmf_res <- nmf(x = m, rank = 3:10, method = "snmf/r", nrun = 10 , .opt = "v", .pbackend = NA)
    
    print("Saving")
    st <- Sys.time()
    saveRDS(object = nmf_res, file = paste0(save_dir, sname, ".RDS"))
    et <- Sys.time()
    print(et - st)
    
    print(paste0("******************************* ", sname, " - end, i = ", i, " *******************************"))
  }
}
end_time <- Sys.time()
end_time - start_time

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Set-up patient pairs
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Using only cells from confidence levels 1 and 2
malignant_cells <- meta_data %>% filter(CellType == "Malignant", ClassConf %in% c("1", "2"))

pt_pairs <- patient_pairs(malignant_cells)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Derive NMF meta-programs for the malignant cells
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

####################################################################################################################################
# Load results
####################################################################################################################################

nmf_res_dir <- paste0(nmf_res_path, "/Malignant/")

nmf_files <- list.files(nmf_res_dir)

Genes_nmf_w_basis <- lapply(nmf_files, function(x) {
  f <- readRDS(file = paste0(nmf_res_dir, x))
  
  res <- get_nmf_programs(f$fit, x)
  
  rm(f)
  gc()
  
  return(res)
})
names(Genes_nmf_w_basis) <- nmf_files

Genes_nmf_w_basis <- Genes_nmf_w_basis[get_sname(names(Genes_nmf_w_basis)) %in% pt_pairs$Sample]

####################################################################################################################################
# Analyze NMF results
####################################################################################################################################

# This will run the NMF meta-program derivation algotirhm
malignant_nmf_metaprograms <- derive_NMF_metaprograms(Genes_nmf_w_basis = Genes_nmf_w_basis, verbose = T, n_genes = 50)

MP <- malignant_nmf_metaprograms$MP_tbl

MP_list <- malignant_nmf_metaprograms$MP_list

# This will run gene-set enrichment analysis on each MP to aid in annotation
mp_en <- metaprograms_enrichment(mp_list = MP_list, pathways = default_pathways(), min_gs_size = 40, max_gs_size = 250)

mp_en %>%
  arrange(MP, desc(MP_Freq)) %>%
  group_by(MP) %>%
  top_n(n = 10, wt = MP_Freq) %>%
  View()

# MP16 includes genes from a mixture of other MPs and was excluded
MP_list <- MP_list[-16]
MP <- MP[, -16]

# Following ispection this is how the MPs were annotated
MP2FUNCTION <- c(MP_1 = "RP", MP_2 = "OPC", MP_3 = "CC", MP_4 = "AC", MP_5 = "Hypoxia",
                 MP_6 = "MES", MP_7 = "NPC", MP_8 = "GPC", MP_9 = "ExN", MP_10 = "Stress1",
                 MP_11 = "MIC", MP_12 = "LQ", MP_13 = "Cilia", MP_14 = "NRGN", MP_15 = "Stress2")

MP_list_named <- setNames(MP_list, paste0(names(MP_list), "_", MP2FUNCTION[names(MP_list)]))

MP_named <- MP
colnames(MP_named) <- paste0(colnames(MP), "_", MP2FUNCTION[colnames(MP)])

malignant_nmf_metaprograms$MP2FUNCTION <- MP2FUNCTION

malignant_nmf_metaprograms$malignant_cells <- malignant_cells

saveRDS(malignant_nmf_metaprograms, paste0(DATA_ROOT, "nmf_metaprograms_malignant_idhwt.RDS"))
write.csv(x = MP_named, file = paste0(TABLES_ROOT, "MP_malignant_named.csv"), col.names = T, row.names = F)

####################################################################################################################################
# Score the cells for the MPs
####################################################################################################################################

# This will score the cells in each sample for all MPs
MP_scores <- score_within_samples(umi_data_list = umi_data_all, md = malignant_cells, sigs = MP_list_named)

saveRDS(MP_scores, paste0(DATA_ROOT, "MP_scores_Malignant_idhwt.RDS"))

MP_scores$MaxMP <- apply(MP_scores[, names(MP_list_named[-length(MP_list_named)])], 1, function(x) names(MP_list_named[-length(MP_list_named)])[which.max(x)])
table(MP_scores$MaxMP)

table(MP_scores$Timepoint, MP_scores$MaxMP)
table(MP_scores$Timepoint, MP_scores$MaxMP) / rowSums(table(MP_scores$Timepoint, MP_scores$MaxMP)) * 100

####################################################################################################################################
# Plot some diagnostics
####################################################################################################################################

sampled_cells <- MP_scores %>% 
  group_by(Patient, Timepoint) %>% 
  sample_n(15) %>%
  ungroup()

sampled_cells <- sampled_cells %>%
  dplyr::select(CellID, Sample, Patient, Timepoint, starts_with("MP"), MaxMP)

scores <- sampled_cells %>% dplyr::select(starts_with("MP"))
scores <- as.matrix(scores)
rownames(scores) <- sampled_cells$CellID
scores <- t(scores)

Heatmap(scores[-16, ], use_raster = F, raster_by_magick = T,
        col = circlize::colorRamp2(c(-.5, -.25, 0, .25, .5), c("blue", "dodgerblue", "white", "orange", "red")),
        show_row_names = T, show_column_names = F,
        cluster_columns = T, cluster_rows = T,
        clustering_distance_rows = "pearson", clustering_method_rows = "complete",
        clustering_distance_columns = "pearson", clustering_method_columns = "complete",
        column_title = "Cells", column_title_gp = gpar(fontsize = 20),
        row_title = "MP scores", row_title_gp = gpar(fontsize = 20),
        top_annotation = HeatmapAnnotation(Timepoint = setNames(sampled_cells$Timepoint, sampled_cells$CellID),
                                           col = list(State = setNames(brewer_pal(palette = "Set1")(length(levels(sampled_cells$State_all))), levels(sampled_cells$State_all)),
                                                      Timepoint = c("T1" = "#66C2A5", "T2" = "#FC8D62", "T3" = "#8DA0CB")),
                                           show_annotation_name = FALSE,
                                           annotation_legend_param = list(title_gp = gpar(fontsize = 20),
                                                                          labels_gp = gpar(fontsize = 14), border = "black",
                                                                          grid_height = unit(7, "mm"))),
        heatmap_legend_param = list(at = c(-.5, -.25, -0, .25, .5), title = "Score",
                                    labels = c("-.5", "-.25", "0", ".25", ".5"),
                                    title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 14), title_position = "topleft", border = "black",
                                    grid_height = unit(7, "mm")),
        column_title_side = "bottom",
        column_split = as.factor(sampled_cells$MaxMP), border = T,
        row_gap = unit(0,"mm"), cluster_column_slices = T, cluster_row_slices = F)

hc <- fastcluster::hclust(as.dist(cor(t(scores[-16, ]))), method = "average")

dm <- melt(scores) %>%
  as_tibble()
colnames(dm) <- c("MP", "CellID", "Score")

dm$MP <- factor(dm$MP, hc$labels[hc$order])

dm <- dm %>%
  left_join(sampled_cells %>%
              select(CellID, MaxMP), by = "CellID")

dm$MaxMP <- gsub("MP_.*_", "", dm$MaxMP)
dm$MaxMP <- factor(dm$MaxMP, gsub("MP_.*_", "", hc$labels[hc$order]))

ggplot(dm, aes(x = CellID, y = MP, fill = Score)) +
  facet_grid(cols = vars(MaxMP), scales = "free") +
  geom_raster() +
  scale_fill_gradient2(name = "", low = "dodgerblue", mid = "white", high = "red",
                       labels = c("-1", "", "0", "", "1"), limits = c(-1, 1), oob = squish) +
  xlab("Cells") +
  ylab("MP scores") +
  theme_gbm_pvsr(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size = 16))

cor_m <- cor(MP_scores %>% select(starts_with("MP"), - MP_16_GlioNeural))

Heatmap(cor_m, use_raster = F, raster_by_magick = T,
        col = circlize::colorRamp2(c(-.5, -.25, 0, .25, .5), c("blue", "dodgerblue", "white", "orange", "red")),
        show_row_names = T, show_column_names = T,
        cluster_columns = T, cluster_rows = T,
        clustering_distance_rows = "pearson", clustering_method_rows = "complete",
        clustering_distance_columns = "pearson", clustering_method_columns = "complete",
        column_title = "MP scores", column_title_gp = gpar(fontsize = 20),
        row_title = "MP scores", row_title_gp = gpar(fontsize = 20),
        heatmap_legend_param = list(at = c(-.5, -.25, -0, .25, .5), title = "Score",
                                    labels = c("-.5", "-.25", "0", ".25", ".5"),
                                    title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 14), title_position = "topleft", border = "black",
                                    grid_height = unit(7, "mm")),
        column_title_side = "bottom",
        row_gap = unit(0,"mm"), cluster_column_slices = F, cluster_row_slices = F)
