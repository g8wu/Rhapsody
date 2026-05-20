# Single Seurat Object
#system("pip install leidenalg")
# install.packages("devtools")
# devtools::install_github("immunogenomics/presto")

# LJI Seurat v5 module adjust
#remotes::install_version("matrixStats", version="1.1.0")
#.rs.restartR()
set.seed(99)
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if(!require("remotes", quietly = T)) install.packages("remotes")
# if (!require("harmony", quietly = TRUE)) BiocManager::install("harmony", force = T)
library(harmony)
library(tidyverse)
library(dplyr)
library(scDblFinder)
library(sctransform)
library(ggplot2)
library(purrr)
library(leidenAlg)
library(Seurat)

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

# Input Settings ############# 
project <- "C104"
file <- ".rds"
alg <- 4
resRange <- seq(0.2, 0.5, by = 0.1)
setwd("/mnt/BioAdHoc/Groups/Collaborators/ben.croker/dress/PBMC")
print("_______________________________________________")

# Global variables
rnaPCs <- 30
adtPCs <- 20
algKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")
print("_______________________________________________")

# For laptops/PCs override default variable size limit
#n <- 50     # Let variables reach up to n GB
#options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB

# Read in files ###########
rds <- readRDS(file)

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
#            "C56-NCAM16.2", "CD62L", "CD63", "CD64-FCGR1A", "CD86", "CD95-FAS", "CXCR2",
#            "Siglec8", "CD66b")

# COCCI
# abseq <- c("CD105", "CD115", "CD117", "CD11a", "CD11b", "CD11c", "CD162", "CD16/32",
#            "CD184", "CD19", "CD31", "CD326", "CD335", "CD41", "CD45", "CD48",
#            "CD62L", "CD71", "CXCR2", "Clec7a", "CD150", "F4/80", "Ly6A/E", "Ly6G",
#            "NK1.1", "SiglecF", "TCRb", "TER119")
# DUPI
abseq <- c("CD101", "CD10", "CD11b", "CD11c", "CD123-IL3RA", "CD14", "CD15-FUT4",
           "CD162-SELPLG", "CD16", "CD183-CXCR3", "CD184-CXCR4", "CD193-CCR3",
           "CD194-CCR4", "CD19", "CD32", "CD33", "CD34", "CD3", "CD41-ITGA2B",
           "CD44", "CD56", "CD62L", "CD63", "CD64", "CD86", "CD95", "CXCR2",
           "Siglec8", "CD66b", "FCER1A", "HLA-DR-CD74")


print(" ----------------------------------------- Abseq markers:")
rownames(rds@assays$ADT)
print(" ----------------------------------------- Rename To:")
print(abseq)

# Apply QC thresholds ############
nFeatLower <- 50
nFeatUpper <- 3000
nCountLower <- 100
nCountUpper <- 8000
mtPct <- 20

# Calc mito pct and Rename ADT
temp <- rds$ADT@counts
rownames(temp) <- abseq_clean
rds[["ADT"]] <- CreateAssayObject(counts = temp)


rds[["percent.mt"]] <- PercentageFeatureSet(rds, pattern = "(?i)^mt-")

# QC Violinplots
rds@project.name <- project
features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
pdf(paste(project, "QC VlnPlot.pdf"))
print(VlnPlot(rds, features = features, ncol = length(features), group.by = "orig.ident"))
dev.off()

# Set QC Thresholds
rds <- subset(rds, subset = nFeature_RNA > nFeatLower & nFeature_RNA < 
                      nFeatUpper & nCount_RNA > nCountLower & nCount_RNA < 
                      nCountUpper & percent.mt < mtPct)
pdf(paste(project, "postQC VlnPlot.pdf"))
print(VlnPlot(rds, features = features, ncol = length(features), group.by = "orig.ident"))
dev.off()

# RNA Normalization ############## 
DefaultAssay(rds) <- 'RNA'
rds <- SCTransform(rds, vars.to.regress = "percent.mt")
rds <- FindVariableFeatures(rds, assay = "RNA", nfeatures = 3000)
rds <- ScaleData(rds)
rds <- RunPCA(rds, npcs = rnaPCs)
rds <- RunUMAP(rds, reduction = "pca", reduction.name = "rna.umap",
                       reduction.key = 'rnaUMAP_', dims = 1:rnaPCs)
print(paste(Sys.time(), "RNA SCT norm done"))

# ADT Normalization ############## 
DefaultAssay(rds) <- 'ADT'
abseq <- rownames(rds)
rdsA <- NormalizeData(rds, normalization.method = 'CLR', margin = 2)
rdsA <- ScaleData(rdsA, features = abseq)
rdsA <- FindVariableFeatures(rdsA, features = rownames(rdsA[["ADT"]]))
rdsA <- RunPCA(rdsA, npcs = adtPCs, reduction.name = 'pca', approx = FALSE)
rdsA <- RunUMAP(rdsA, reduction = "pca", reduction.name = "adt.umap",
                       reduction.key = 'adtUMAP_', dims = 1:adtPCs)

print(paste(Sys.time(), "ADT Norm done"))

# INTEGRATION ############## 
rds[["ADT"]] <- rdsA[["ADT"]]
rds[["apca"]] <- rdsA[["pca"]]
rds[["adt.umap"]] <- rdsA[["adt.umap"]]

DefaultAssay(rds) <- "SCT"

# scDoubletFinder ####
sce <- as.SingleCellExperiment(rds, assay = "SCT")
sce <- scDblFinder(sce, samples = "orig.ident")

rds$scDblFinder.class <- sce$scDblFinder.class
rds$scDblFinder.class <- sce$scDblFinder.score

rds <- AddMetaData(rds, metadata=sce$scDblFinder.score, col.name='scDblFinder_score')
rds <- AddMetaData(rds, metadata=sce$scDblFinder.class, col.name='scDblFinder_class')

pdf(paste0(rds@project.name, "-UMAPRNA_scDblFinder.pdf"))
DimPlot(rds, reduction="rna.umap", group.by="scDblFinder_class", cols=c('grey', 'red'), order=TRUE)
FeaturePlot(rds, reduction="rna.umap", features='scDblFinder_score', order=TRUE)
dev.off()

prop.table(table(rds$scDblFinder_class))*100
summary(rds$nCount_RNA[rds$scDblFinder_class=='doublet'])
summary(rds$nCount_RNA[rds$scDblFinder_class=='singlet'])

#saveRDS(rds, paste0(rds@project.name, "DbFinder.rds"))
Idents(rds) <- rds$scDblFinder_class
rds <- subset(rds, idents= "singlet")

##############  RNA ONLY WNN ############## 
## WNN data: combinedRNA[['weighted.nn']]
## WNN graph: combinedRNA[["wknn"]],
## SNN graph used for clustering: combinedRNA[["wsnn"]]
## Cell-specific modality weights: combinedRNA$RNA.weight
print("_______________________________________________ Starting RNA WNN Clustering")
rds <- FindVariableFeatures(rds, assay = "SCT", selection.method = "vst", nfeatures = 3000)
rds <- FindNeighbors(rds, dims = 1:30, graph.name = "RNA_nn")
rds <- RunUMAP(rds, dims = 1:30)

# Save metadata
rds@project.name <- paste0(project, "-RNA")
rds@misc <- list(algKey[alg], "umap", "RNA")

# Multi-res Clustering
pdf(paste0(rds@project.name, "-RNAWNNUMAP.pdf"))
for (feat in features){
  print(FeaturePlot(rds, features = feat, raster = F))
for(res in resRange) {
  print(res)
  rds <- FindClusters(rds, algorithm = alg, resolution = res, graph.name = "RNA_nn")
  Idents(rds) <- rds[[paste0("RNA_nn_res.", res)]][,1]
  print(DimPlot(rds, reduction = 'rna.umap', label = TRUE, raster = F) +
          labs(title = paste0("RNA WNN UMAP ", algKey[alg],": ", res)))
  }
}
dev.off()
saveRDS(rds, file= paste0(rds@project.name, "-RNAClust.rds"))

print(paste0(Sys.time(), " -> RNA ONLY multi res WNN done!"))

