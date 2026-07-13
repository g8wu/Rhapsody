set.seed(99)
#.rs.restartR()
#system("pip install leidenalg")
install.packages("devtools")
devtools::install_github("immunogenomics/presto")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#if (!require("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!require("BiocManager", quietly = TRUE)) BiocManager::install("glmGamPoi")
if(!require("remotes", quietly = T)) install.packages("remotes")
if (!require("harmony", quietly = TRUE)) BiocManager::install("harmony")
if (!require("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")
if (!require("ADTnorm", quietly = TRUE)) remotes::install_github("yezhengSTAT/ADTnorm", build_vignettes = FALSE)
library(packrat)
library(harmony)
library(tidyverse)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(sctransform)
library(ggplot2)
library(purrr)
library(reticulate)
# Install python and make virtual env with pythom -m venv VenvName
# use_virtualenv("~/VENV")
# py_config()
library(leidenAlg)
library(ADTnorm) #2
library(Seurat)
options(Seurat.object.assay.version = 'v5')

# Bash Script Settings ############# 
# FindCluster alg 1:  Louvain | fast & effective but not for complex datasets,
#                 2:  Refined Louvain | multilevel refined clusters but computationally heavier,
#                 3:  Smart Local Moving (SLM) | Louvain w/ more granularity, slower but good for complex datasets
#                 4:  Leiden | fast & more accurate w/well-connected clusters, requires leidenalg Python package)
# project <- args[1]
# fileType <- args[2]
# alg <- args[3]
# resRange <- seq(args[4], args[5], by = 0.1)
# multRate <- args[6]            # calculated multiplet rate from upstream files
# 
# args <- commandArgs(trailingOnly = TRUE)
# print(paste("Project: ", project))
# print(paste("fileType: ", fileType))
# print(paste("Cluster alg: ", alg))
# print(paste("resRange: ", args[4], args[5]))
# print(paste("MultRate: ", multRate))

# Read from tsv files (Depreciated) ####
# expression_matrix <- ReadMtx(
#   mtx = "GSE291519_matrix.mtx.gz", features = "GSE291519_features.tsv.gz",
#   cells = "GSE291519_barcodes.tsv.gz"
# )
# rds <- CreateSeuratObject(counts = expression_matrix)

# Input Settings ############# 
#sink("log_file.txt")
project <- "Nomid"
fileType <- ".rds"
alg <- 4
resRange <- seq(0.4, 0.6, by = 0.1)
setwd("/mnt/bioadhoc/Groups/Collaborators/ben.croker/nomid")
print("_______________________________________________")

# Global variables
rnaPCs <- 30
adtPCs <- 10
algKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")
print("_______________________________________________")

# For laptops/PCs override default variable size limit
n <- 50     # Let variables reach up to n GB
options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB

# Read in files ###########
wkdir <- getwd()
setwd(wkdir)
rds_files <- list.files(path = paste0(wkdir, "/samples/og"), pattern = fileType, full.names = TRUE)
print(rds_files)
list <- list()
list <- map(rds_files, readRDS)

# Truncate Abseq ###########
# NOMID
# abseq <- c("CD105", "CD115", "CD117", "CD11a", "CD11b", "CD11c", "CD162", "CD16/32",
#            "CD184", "CD19", "CD335", "CD41", "CD45", "CD48", "CD62L", "CD71",
#            "CXCR2", "Clec7a", "CD150", "F4/80", "Ly6A/E", "Ly6G", "NK1.1", "SiglecF",
#            "TCRb", "TER119")
# DRESS
# abseq <- c("CD101", "CD10", "CD11b", "CD11c", "CD123", "CD14", "CD15", "CD162-SELPLG",
#            "CD16-FCGR3A", "CD183-CXCR3", "CD184-CXCR4", "CD193-CCR3", "CD194-CCR4",
#            "CD19", "CD32-FGCR2A", "CD33", "CD34", "CD3", "CD41-ITGA2B", "CD44",
#            "C56-NCAM16.2", "CD62L-SELL", "CD63", "CD64-FCGR1A", "CD86", "CD95-FAS", "CXCR2",
#            "Siglec8", "CD66b", "FCER1A", "HLADR-CD74")

# COCCI
# abseq <- c("CD105", "CD115", "CD117", "CD11a", "CD11b", "CD11c", "CD162", "CD16/32",
#            "CD184", "CD19", "CD31", "CD326", "CD335", "CD41", "CD45", "CD48",
#            "CD62L", "CD71", "CXCR2", "Clec7a", "CD150", "F4/80", "Ly6A/E", "Ly6G",
#            "NK1.1", "SiglecF", "TCRb", "TER119")
# DUPI
# abseq <- c("CD101", "CD10", "CD11b", "CD11c", "CD123-IL3RA", "CD14", "CD15-FUT4", 
#            "CD162-SELPLG", "CD16", "CD183-CXCR3", "CD184-CXCR4", "CD193-CCR3", 
#            "CD194-CCR4", "CD19", "CD32", "CD33", "CD34", "CD3", "CD41-ITGA2B", 
#            "CD44", "CD56", "CD62L", "CD63", "CD64", "CD86", "CD95", "CXCR2", 
#            "Siglec8", "CD66b", "FCER1A", "HLA-DR-CD74")


print(" ----------------------------------------- Abseq markers:")
print(rownames(list[[1]]@assays$ADT))
print(" ----------------------------------------- Rename To:")
print(abseq)

# Apply QC thresholds ############
nFeatLower <- 50
nFeatUpper <- 5000
nCountLower <- 100
nCountUpper <- 18000
mtPct <- 25

features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "ADT_total")
# Calc mito pct and Rename ADT
for (i in seq_along(list)){
  # Trunc abseq names
  temp <- list[[i]]@assays$ADT@counts
  rownames(temp) <- abseq
  list[[i]][["ADT"]] <- CreateAssayObject(counts = temp)
  
  # Calc percent mito reads
  list[[i]][["percent.mt"]] <- PercentageFeatureSet(list[[i]], pattern = "(?i)^mt-")
  
  # Calc total abseq reads
  list[[i]]$ADT_total <- colSums(list[[i]]@assays$ADT@data)
  
  # Apply QC thresholds
  list[[i]] <- subset(list[[i]], subset = nFeature_RNA > nFeatLower & nFeature_RNA < 
                        nFeatUpper & nCount_RNA > nCountLower & nCount_RNA < 
                        nCountUpper & percent.mt < mtPct)
}
rm(temp)

#PRINT
pdf(paste(project, "VlnPlotQC.pdf", sep = "-"))
for (i in seq_along(features)){
  # Initialize variable to store all data
  plot_data <- data.frame()
  for (j in seq_along(list)){
    # extract QC data from each list object
    metaData <- data.frame(Sample = list[[j]]@meta.data$orig.ident,
                           Feat = list[[j]]@meta.data[[features[i]]])
    # Append to full dataset
    plot_data <- bind_rows(plot_data, metaData)
  }
  print(ggplot(plot_data, aes(x = Sample, y = Feat, fill = Sample)) +
          geom_violin() + theme_minimal() + 
          theme(axis.text.x = element_text(angle = 90)) +
          labs(title = features[i], x = "Sample", y = features[i]))
}
dev.off()

# RNA Normalization ############## 
print("____________________________________________ RNA SCT nomalization start")
list <- lapply(list, function(i) {
  DefaultAssay(i) <- 'RNA'
  i <- SCTransform(i, vars.to.regress = "percent.mt")
  return(i)
})
list
print(paste(Sys.time(), "RNA SCT norm done"))
gc()

# RNA Integration ############## 
print("____________________________________________ RNA integration start")
features <- SelectIntegrationFeatures(list, nfeatures = 3000)
list <- PrepSCTIntegration(list, anchor.features = features)
list <- lapply(list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(list, anchor.features = features, dims = 1:rnaPCs,
                                  normalization.method = "SCT")
combinedRNA <- IntegrateData(anchors, normalization.method = "SCT")
DefaultAssay(combinedRNA) <- "integrated"
gc()

combinedRNA <- RunPCA(combinedRNA, npcs = rnaPCs, reduction.name = "rpca")

## RNA BeC ####
combinedRNA <- RunHarmony(combinedRNA, group.by.vars = "orig.ident", 
                          reduction.use = "rpca", reduction.save = "harmony.rpca")
# With BeC
combinedRNA <- RunUMAP(combinedRNA, reduction = "harmony.rpca", reduction.name = "harmony.rumap",
                       dims = 1:rnaPCs)
# Without BeC
combinedRNA <- RunUMAP(combinedRNA, reduction = "rpca", reduction.name = "rumap",
                       dims = 1:rnaPCs)

print(paste0(Sys.time(), " -> RNA SCTransform and integration done and saved!"), append = TRUE)

# ADT Normalization ############## 
# Save CLR normalized and merged to separate object (preserves compositional nature but no batch effect correction)
# Save logNormalized and integrated for comparison?
print("____________________________________________ ADT CLR Normalization")
abseq <- rownames(list[[1]]@assays$ADT)
list <- lapply(list, function(i) {
  DefaultAssay(i) <- 'ADT'
  i <- NormalizeData(i, normalization.method = 'CLR', margin = 2)
  i <- ScaleData(i)
  i <- RunPCA(i, reduction.name = "apca")
  return(i)
})
combinedADT <- merge(list[[1]], y = list[-1])
DefaultAssay(combinedADT) <- "ADT"

## ADTnorm ####
# Get raw ADT counts with cell ID rows, Abseq features cols
cell_x_adt <- t(combinedADT@assays$ADT@counts)
cell_x_feature <- data.frame(combinedADT@active.ident)

# Make sure there is a 'sample' column set to an ident
colnames(cell_x_feature) <- 'sample'
head(cell_x_feature)

cell_x_adt_norm = ADTnorm(
  cell_x_adt = cell_x_adt, 
  cell_x_feature = cell_x_feature, 
  save_outpath = getwd(), 
  study_name = combinedADT@project.name,
  marker_to_process = newAbseq,
  bimodal_marker = NULL,             # default NULL: try different settings to find bimodal peaks
  trimodal_marker = c("CD45"),       # CD4 and CD45RA tend to have 3 peaks
  # setting the CD3 uni-peak of buus_2021_T study to positive peak if only one peak is detected for CD3 marker
  # positive_peak = list(ADT = "CD3", sample = "buus_2021_T"), 
  positive_peak = list(ADT = "CD19"), 
  brewer_palettes = "Dark2",
  save_fig = TRUE,
  target_landmark_location = "fixed",
  shoulder_valley = T,               # Look for "shoulder" as pos peak (technical variation -> no clear separation b/w neg/pos)
  #multi_sample_per_batch = T,        # Omit aligning the one pos peak
  #customize_landmark = T             # Manual adjustment UI 
)

# Put ADTNorm matrix back into rds
# combinedADT@assays$ADT@data <- t(cell_x_adt_norm)

# Save
# saveRDS(combinedADT, paste0(project, "ADTnorm.rds"))
# gc()
# print(paste(Sys.time(), "ADT Norm done"))

# Only features are from abseq panel
VariableFeatures(combinedADT) <- rownames(combinedADT[['ADT']])
combinedADT <- ScaleData(combinedADT)
combinedADT <- RunPCA(combinedADT, npcs = length(rownames(rds@assays$ADT)), reduction.name = 'apca')
combinedADT <- RunUMAP(combinedADT, reduction = "apca", reduction.name = "aumap", 
                       dims = 1:adtPCs)
print(paste0(Sys.time(), " -> ADT & RNA independent normalization done!"))

# INTEGRATION ############## 
print("____________________________________________ INTEGRATION")
combinedRNA[["ADT"]] <- combinedADT[["ADT"]]
combinedRNA[["apca"]] <- combinedADT[["apca"]]
combinedRNA[["aumap"]] <- combinedADT[["aumap"]]

#PRINT
pdf(paste(project, "harmonyPCAUMAP.pdf", sep = "-"))
ElbowPlot(combinedRNA, reduction = "rpca")
DimPlot(combinedRNA, reduction = "rpca", group.by = "orig.ident")
DimPlot(combinedRNA, reduction = "rumap", group.by = "orig.ident")
DimPlot(combinedRNA, reduction = "harmony.rpca", group.by = "orig.ident")
DimPlot(combinedRNA, reduction = "harmony.rumap", group.by = "orig.ident")
ElbowPlot(combinedRNA, reduction = "apca")
DimPlot(combinedRNA, reduction = "apca", group.by = "orig.ident")
DimPlot(combinedRNA, reduction = "aumap", group.by = "orig.ident")
dev.off()

# SAVE!
combinedRNA@project.name <- project
saveRDS(combinedRNA, file = paste(project, "integrated.RDS", sep = "-"))
print(paste0(Sys.time(), " -> RNA & ADT INTEGRATION successful, saved before WNN!"))
DefaultAssay(combinedRNA) <- "integrated"
rds <- combinedRNA
rm(combinedADT)
rm(combinedRNA)
gc()

# scDoubletFinder ####
sce <- as.SingleCellExperiment(rds, assay = "SCT")
sce <- scDblFinder(sce, samples = "orig.ident")

rds$scDblFinder.class <- sce$scDblFinder.class
rds$scDblFinder.class <- sce$scDblFinder.score

rds <- AddMetaData(rds, metadata=sce$scDblFinder.score, col.name='scDblFinder_score')
rds <- AddMetaData(rds, metadata=sce$scDblFinder.class, col.name='scDblFinder_class')

pdf(paste0(project, "-RNA_dbFinderUMAP.pdf"))
DimPlot(rds, reduction="rna.umap", raster = F, group.by="scDblFinder_class", cols=c('grey', 'red'), order=TRUE)
FeaturePlot(rds, reduction="rna.umap", features='scDblFinder_score', raster = F, order=TRUE)
dev.off()

# Save summary
write.csv(data.frame(prop.table(table(rds$scDblFinder_class))*100), paste0(project, "-DblFinder.csv"), row.names = F)
write.table(data.frame(Class = "Doublet", 
                       Statistic = names(summary(rds$nCount_RNA[rds$scDblFinder_class=='doublet'])),
                       Value = as.numeric(summary(rds$nCount_RNA[rds$scDblFinder_class=='doublet']))),
            paste0(project, "-DblFinder.csv"), sep = ",", row.names = F, col.names = F, append = T)
write.table(data.frame(Class = "Singlet", 
                       Statistic = names(summary(rds$nCount_RNA[rds$scDblFinder_class=='singlet'])),
                       Value = as.numeric(summary(rds$nCount_RNA[rds$scDblFinder_class=='singlet']))),
            paste0(project, "-DblFinder.csv"), sep = ",", row.names = F, col.names = F, append = T)
# Save
saveRDS(rds, paste0(rds@project.name, "-DbFind.rds"))
Idents(rds) <- rds$scDblFinder_class
rds <- subset(rds, idents= "singlet")

###########  RNA ONLY  ############## 
## WNN data: combinedRNA[['weighted.nn']]
## WNN graph: combinedRNA[["wknn"]],
## SNN graph used for clustering: combinedRNA[["wsnn"]]
## Cell-specific modality weights: combinedRNA$RNA.weight
print("_______________________________________________ Starting RNA WNN Clustering")
#rds <- FindVariableFeatures(rds, assay = "SCT", selection.method = "vst", nfeatures = 3000)
rds <- FindNeighbors(rds, dims = 1:rnaPCs, reduction = "harmony.rpca")
rds <- RunUMAP(rds, dims = 1:rnaPCs, reduction.name = "harmony.rumap", reduction = "harmony.rpca")

# Save metadata
rds@misc$clustAlg <- algKey[alg]
rds@misc$umap <- "harmony.rumap"

# Multi-res Clustering
pdf(paste0(rds@project.name, "-UMAPRNA.pdf"))
for(res in resRange) {
  print(res)
  rds <- FindClusters(rds, algorithm = alg, resolution = res, graph.name = "integrated_snn", random.seed = 1)
  print(DimPlot(rds, reduction = 'harmony.rumap', label = TRUE, raster = F) +
          labs(title = paste0("RNA WNN UMAP ", algKey[alg],": ", res)))
}
dev.off()
# saveRDS(rds, file= paste0(rds@project.name, "-RNAClust.rds"))

print(paste0(Sys.time(), " -> RNA ONLY multi res WNN done!"))

############RNA & ADT ############## 
# Use RNA UMAP and ADT PCA for MultiModal
print("_______________________________________________ Starting ADT & RNA WNN Clustering")
rds <- FindMultiModalNeighbors(rds, reduction.list = list("harmony.rpca", "apca"),
                               dims.list = list(1:rnaPCs, 1:adtPCs))
rds <- RunUMAP(rds, nn.name = "weighted.nn", reduction.name = "wnn.umap")
# Save metadata 
rds@misc$clustAlg <- algKey[alg]

# Multi-res Clustering
pdf(paste0(rds@project.name, "-UMAPWNN.pdf"))
for(res in resRange) {
  print(res)
  rds <- FindClusters(rds, algorithm = alg, resolution = res, graph.name = "wsnn")
  print(DimPlot(rds, reduction= "wnn.umap", label = T, raster = F) + 
          labs(title = paste0("ADT&RNA UMAP ", algKey[alg],": ", res)))
}
dev.off()
saveRDS(rds, file= paste0(rds@project.name, "-RNACLR.rds"))

print(paste0(Sys.time(), " -> RNA & ADT multi res WNN done!"))

