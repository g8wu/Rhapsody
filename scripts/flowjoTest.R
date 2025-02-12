########################################
# Run 1 or few samples, includes FloJo conversion
# 3/25/2024 Gio Wu
# Refs:
# Based off Basic Protocols 2024, Li et al.
# https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
########################################
set.seed(99)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
library(harmony)
library(gridExtra)
library(rstudioapi)
library(patchwork)
library(S4Vectors)
library(tidyverse)
library(sctransform)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(flowCore)
library(Biobase)
library(rstudioapi)
library(flowCore)

project <- "314"
wkdir <- paste0(getwd(), "/Projects/")
pre <- readRDS(paste0(wkdir, project, "/pre", project, ".rds"))
post <- readRDS(paste0(wkdir, project, "/post", project, ".rds"))
list <- list(pre, post)

# QC
prePCT_MT <- 20
preNFEAT_UP <- 3000
preNFEAT_LOW <- 200
preNCOUNT_UP <- 30000
preNUM_PCA_DIM <- 30
preNUM_ADT_PCA_DIM <- 18

postPCT_MT <- 18
postNFEAT_UP <- 4400
postNFEAT_LOW <- 200
postNCOUNT_UP <- 40000
postNUM_PCA_DIM <- 30
postNUM_ADT_PCA_DIM <- 18

features = c("nFeature_RNA", "nCount_RNA", "percent.mt") 
pre[["percent.mt"]] <- PercentageFeatureSet(pre, pattern = "^MT-") 
post[["percent.mt"]] <- PercentageFeatureSet(post, pattern = "^MT-") 
dim(pre)
dim(post)
VlnPlot(pre, features = features, pt.size = 0.01, ncol = 3) & labs(x="")
VlnPlot(post, features = features, pt.size = 0.01, ncol = 3) & labs(x="")

# SET QC PARAMS
pre <- subset(pre, subset = nFeature_RNA > preNFEAT_LOW & nFeature_RNA < preNFEAT_UP & percent.mt < prePCT_MT & nCount_RNA < preNCOUNT_UP)
post <- subset(post, subset = nFeature_RNA > postNFEAT_LOW & nFeature_RNA < postNFEAT_UP & percent.mt < postPCT_MT & nCount_RNA < postNCOUNT_UP)
dim(pre)
dim(post)
VlnPlot(pre, features = features, pt.size = 0.01, ncol = 3) & labs(x="")
VlnPlot(post, features = features, pt.size = 0.01, ncol = 3) & labs(x="")

saveRDS(pre, file = paste0(wkdir, project, "/pre", project, "_QC.rds"))
cat(sprintf("\n ==>> QC file saved for pre"))
saveRDS(pre, file = paste0(wkdir, project, "/post", project, "_QC.rds"))
cat(sprintf("\n ==>> QC file saved for post"))

# FLOWFRAME FCS file
data_matrix <- GetAssayData(pre, assay = "WTA", layer = "data")
data_matrix <- as.matrix(GetAssayData(pre, assay = "WTA", slot = "data", layer = "data"))
ff <- flowFrame(exprs = data_matrix)
# Save the flowFrame as an FCS file
write.FCS(ff, paste0(wkdir, project, "/pre", project, ".FCS"))
cat(sprintf("\n ==>> preFCS file saved"))

data_matrix <- GetAssayData(post, assay = "WTA", layer = "data")
data_matrix <- as.matrix(GetAssayData(post, assay = "WTA", slot = "data", layer = "data"))
ff <- flowFrame(exprs = data_matrix)
# Save the flowFrame as an FCS file
write.FCS(ff, paste0(wkdir, project, "/post", project, ".FCS"))
cat(sprintf("\n ==>> postFCS file saved"))

save(fcs_demo_qc, file = paste0(wkdir, project, "/post", project, ".FCS"))
meta_list <- c("nCount_RNA", 
               "nFeature_RNA", 
               "percent.mt", 
               "seurat_clusters", 
               "m3_pseudotime")

func_create_fcs(input_seurat = pre, include_protein_assay_name = "ABSEQ", 
                meta_data_to_pull = 
                outfile_name = paste0("pre", project, ".FCS"))
