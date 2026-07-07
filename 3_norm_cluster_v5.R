# Libraries ####
set.seed(99)
#.rs.restartR()
# system("pip install leidenalg")
# install.packages("devtools")
# devtools::install_github("immunogenomics/presto")
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("glmGamPoi")
# if(!require("remotes", quietly = T)) install.packages("remotes")
# if (!require("harmony", quietly = TRUE)) BiocManager::install("harmony")
# if(!require("scCustomize", quietly = T)) install.packages("scCustomize")
# if (!require("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")
# if (!require("ADTnorm", quietly = TRUE)) { BiocManager::install(c("cytolib","flowCore", "ncdfFlow", "flowViz", "flowWorkspace","flowStats"))
#   remotes::install_github("yezhengSTAT/ADTnorm", build_vignettes = FALSE)
# }
# if (!require("leidenAlg", quietly = T)) {install.packages("leidenAlg")}
library(harmony)
library(tidyverse)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(sctransform)
library(ggplot2)
library(purrr)
library(leidenAlg)
library(ADTnorm)
library(Seurat)

# Input Settings ####
#sink("log_file.txt")
rds <- readRDS("RIME-postQC.rds")
alg <- 4
resRange <- seq(0.1, 0.5, by = 0.1)
# setwd("/mnt/bioadhoc/Groups/Collaborators/ben.croker/nomid")
print("_______________________________________________")

# Global variables
rnaPCs <- 50
adtPCs <- 20
algKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")

# For laptops/PCs override default variable size limit
n <- 50     # Let variables reach up to n GB
options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB

# RNA Normalization ############## 
print("____________________________________________ SCTransform start")
DefaultAssay(rds) <- 'RNA'
list <- SplitObject(rds, split.by = "orig.ident")
list <- lapply(list, \(x){
  SCTransform(x, assay = 'RNA', 
              new.assay.name = "SCT",
              vars.to.regress = "percent.mt")
})
print(paste(Sys.time(), "RNA SCT norm done"))
gc()

# SAVE!!!
saveRDS(rds, paste0(rds@project.name, "-SCT.rds"))

# ADTnorm #####
print("____________________________________________ ADTnorm start")
cell_x_adt <- t(rds@assays$ADT@counts)
cell_x_feature <- data.frame(rds@active.ident)

# Make sure there is a 'sample' column set to an ident
colnames(cell_x_feature) <- 'sample'
head(cell_x_feature)

cell_x_adt_norm = ADTnorm(
  cell_x_adt = cell_x_adt, 
  cell_x_feature = cell_x_feature, 
  save_outpath = getwd(), 
  study_name = rds@project.name,
  marker_to_process = rownames(rds@assays$ADT),
  bimodal_marker = NULL,             # default NULL: try different settings to find bimodal peaks
  # trimodal_marker = c("CD45-F11-Ptprc"),       # CD4 and CD45RA tend to have 3 peaks
  # setting the CD3 uni-peak of buus_2021_T study to positive peak if only one peak is detected for CD3 marker
  # positive_peak = list(ADT = "CD3", sample = "buus_2021_T"), 
  positive_peak = list(ADT = "CD3"), 
  brewer_palettes = "Dark2",
  save_fig = TRUE,
  target_landmark_location = "fixed",
  shoulder_valley = T,               # Look for "shoulder" as pos peak (technical variation -> no clear separation b/w neg/pos)
  #multi_sample_per_batch = T,        # Omit aligning the one pos peak
  #customize_landmark = T             # Manual adjustment UI 
)

# Put ADTNorm matrix back into rds
rds@assays$ADT@scale.data <- t(cell_x_adt_norm)
print(paste(Sys.time(), "ADTnorm done"))


# PCA ####
rds <- RunPCA(rds, assay = "SCT", reduction.name = "pca_rna")
VariableFeatures(rds@assays$ADT) <- rownames(rds@assays$ADT)
rds <- RunPCA(rds, assay = "ADT", reduction.name = "pca_adt", npcs = adtPCs)

# Harmonize ####
# Batch effect correct (BeC)
rds <- RunHarmony(rds, group.by.vars = "orig.ident", 
                  reduction = "pca_rna", 
                  reduction.save = "harmony_rna")
rds <- RunHarmony(rds, group.by.vars = "orig.ident", 
                  reduction = "pca_adt", 
                  reduction.save = "harmony_adt")
# UMAP RNA ####
# UMAP Without Harmony
rds <- RunUMAP(rds, reduction = "pca_rna", dims = 1:rnaPCs, reduction.name = "umap")
# UMAP With Harmony
rds <- RunUMAP(rds, reduction = "harmony_rna", dims = 1:rnaPCs, reduction.name = "harmony_rna")
# UMAP With Harmony
rds <- RunUMAP(rds, reduction = "harmony_adt", dims = 1:adtPCs, reduction.name = "harmony_adt")


# Print
pdf(paste0(rds@project.name,"-PCA-Harmony.pdf"))
print(DimPlot(rds, reduction = "umap"))
print(DimPlot(rds, reduction = "harmony_rna"))
print(DimPlot(rds, reduction = "harmony_adt"))
print(DimPlot(rds, reduction = "pca_rna"))
print(DimPlot(rds, reduction = "pca_adt"))
dev.off()

# UMAP RNA + ADT ####
rds <- FindMultiModalNeighbors(rds, reduction.list = list("harmony_rna", "harmony_adt"),
                               dims.list = list(1:rnaPCs, 1:adtPCs))
rds <- RunUMAP(rds, nn.name= "weighted.nn", reduction.name = "wnn.umap")

# DoubletFinder ####
sce <- as.SingleCellExperiment(rds, assay = "SCT")
sce <- scDblFinder(sce, samples = "orig.ident")

rds$scDblFinder.class <- sce$scDblFinder.class
rds$scDblFinder.class <- sce$scDblFinder.score

rds <- AddMetaData(rds, metadata=sce$scDblFinder.score, col.name='scDblFinder_score')
rds <- AddMetaData(rds, metadata=sce$scDblFinder.class, col.name='scDblFinder_class')

pdf(paste0(rds@project.name, "-DbF-UMAP.pdf"))
DimPlot(rds, reduction="harmony_rna.umap", group.by="scDblFinder_class", cols=c('grey', 'red'), order=TRUE)
FeaturePlot(rds, reduction="harmony_rna.umap", features='scDblFinder_score', raster = F, order=TRUE)
DimPlot(rds, reduction="wnn.umap", group.by="scDblFinder_class", cols=c('grey', 'red'), order=TRUE)
FeaturePlot(rds, reduction="wnn.umap", features='scDblFinder_score', raster = F, order=TRUE)
dev.off()

# Save summary
write.csv(data.frame(prop.table(table(rds$scDblFinder_class))*100), paste0(rds@project.name, "-DbF.csv"), row.names = F)
write.table(data.frame(Class = "Doublet", 
                       Statistic = names(summary(rds$nCount_RNA[rds$scDblFinder_class=='doublet'])),
                       Value = as.numeric(summary(rds$nCount_RNA[rds$scDblFinder_class=='doublet']))),
            paste0(rds@project.name, "-DblFinder.csv"), sep = ",", row.names = F, col.names = F, append = T)
write.table(data.frame(Class = "Singlet", 
                       Statistic = names(summary(rds$nCount_RNA[rds$scDblFinder_class=='singlet'])),
                       Value = as.numeric(summary(rds$nCount_RNA[rds$scDblFinder_class=='singlet']))),
            paste0(rds@project.name, "-DblFinder.csv"), sep = ",", row.names = F, col.names = F, append = T)
# Save
saveRDS(rds, paste0(rds@project.name, "Norm-BeC-DbF.rds"))
Idents(rds) <- rds$scDblFinder_class
rds <- subset(rds, idents= "singlet")


# Multi-res Clust ####
## RNA ONLY ####
print("_____________________ Starting RNA ONLY Clustering")
set.seed(99)

rds <- FindNeighbors(rds, dims = 1:rnaPCs, reduction = "harmony_rna")

# Multi-res Clustering
pdf(paste0(rds@project.name, "-UMAP-RNA.pdf"))
for(res in resRange) {
  print(res)
  rds <- FindClusters(rds, algorithm = alg, resolution = res)
  print(DimPlot(rds, reduction = 'harmony_rna.umap', label = TRUE, raster = F) +
          labs(title = paste0("RNA UMAP ", algKey[alg],": ", res)))
}
dev.off()

print(paste0(Sys.time(), " -> RNA ONLY clustering done!"))

## RNA + ADT ####
print("_____________________ Starting RNA & ADT Clustering")
set.seed(99)

# Multi-res Clustering
pdf(paste0(rds@project.name, "-UMAP-RNAADT.pdf"))
for(res in resRange) {
  print(res)
  rds <- FindClusters(rds, algorithm = alg, resolution = res, graph.name = "wsnn")
  print(DimPlot(rds, reduction= "wnn.umap", label = T, raster = F) + 
          labs(title = paste0("ADT&RNA UMAP ", algKey[alg],": ", res)))
}
dev.off()

# SAVE!! ####
saveRDS(rds, file= paste0(rds@project.name, "-RNAadtNorm.rds"))

print(paste0(Sys.time(), " -> RNA & ADT clustering done!"))

# Post Cluster ####
features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "ADT_total")

## QC FeatPlot ####
pdf(paste0(rds@project.name, "-QCUMAP.pdf"))
print(FeaturePlot(rds, reduction = "harmony_umap", features = features, ncol = 2))
print(FeaturePlot(rds, reduction = "wnn.umap", features = features, ncol = 2))
dev.off()

## Select UMAP/res ####
rds@misc$umap <- "harmony_umap"
rds$seurat_clusters <- rds$SCT_snn_res.0.2
Idents(rds) <- "seurat_clusters"

# Continue to 4_dotplot_volcano.R