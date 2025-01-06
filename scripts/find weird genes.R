# find weird genes
set.seed(99)
library(GEOquery)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(ggpubr)
library(Seurat)
library(RColorBrewer)
library(ggrepel)
library(harmony)
library(gridExtra)
library(rstudioapi)
library(patchwork)
library(S4Vectors)
library(tidyverse)
library(sctransform)

# Parameters
pct_MT <- 18
ncount_AB_lower <- 150
ncount_AB_upper <- 100000
nFeature_RNA_lower <- 300
nFeature_RNA_upper <- 5000
nCount_RNA_upper <- 20000
NUM_PCA_DIM <- 30 
NUM_ADT_PCA_DIM <- 18

cart <- "106"
#setwd("C:/Users/Innoscan/Desktop/Gio Wu OneDrive/OneDrive - University of California, San Diego Health/Rhapsody") # Work PC
setwd("C:/Users/Ginny/OneDrive - University of California, San Diego Health/Rhapsody")   # Home PC
savePath <- paste0(cart, "/")

wb_pre <- readRDS(paste0(savePath, "C8-exact_Seurat.rds"))
wb_post <- readRDS(paste0(savePath, "C32-exact_Seurat.rds"))
wb_post <- readRDS(paste0(savePath, "C44_Seurat.rds"))


rds.list <- list(wb_pre, wb_post)

# Split WTA and Abseq, output normalization plots
rds.list <- lapply(X = rds.list, FUN = function(x){
  abseq <- x@assays$RNA$counts[grep("pAbO", rownames(x)),]
  wta <- x@assays$RNA$counts[-(grep("pAbO", rownames(x))),]
  x[["AB"]] <- CreateAssayObject(abseq)
  x[["RNA"]] <- CreateAssayObject(wta)
  DefaultAssay(x) <- "RNA"
  rm(abseq)
  rm(wta)
  return(x)
})

## QC
rds.list <- lapply(X = rds.list, FUN = function(x){
  cat("\n", x$orig.ident[1], "number of cells before filtering : ", ncol(x))
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
  pdf(paste0(savePath, x$orig.ident[1], "_unfiltered_Vlnplot.pdf"))
  print(VlnPlot(x, features = features, pt.size = 0.01, ncol = 3) & labs(x=x$orig.ident[1]))
  testX <- subset(x, subset = nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < pct_MT & nCount_RNA < nCount_RNA_upper)
  dev.off()
  cat("\n", testX$orig.ident[1], "number of cells after filtering : ", ncol(testX))
  pdf(paste0(savePath, x$orig.ident[1], "_filtered_Vlnplot.pdf"))
  print(VlnPlot(testX, features = features, pt.size = 0.01, ncol = 3) & labs(x=x$orig.ident[1]))
  dev.off()
  return(x)
})
cat(sprintf("\n ==>> QC done \n"))

## ADJUST PARAMETERS THEN SET BY RUNNING BELOW
rds.list <- lapply(rds.list, function(x){
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
  x <- subset(x, subset = nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < pct_MT & nCount_RNA < nCount_RNA_upper)
  pdf(paste0(savePath, x$orig.ident[1], "_FINAL_filtered_Vlnplot.pdf"))
  print(VlnPlot(x, features = features, pt.size = 0.01, ncol = 3) & labs(x=""))
  return(x)
})
cat(sprintf("\n ==>> Set QC done \n"))

############ Normalize WTA separately
rds.list <- lapply(rds.list, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
cat(sprintf("\n ==>> Normalization done \n"))

############ Integrate RNA
features <- SelectIntegrationFeatures(rds.list)
rds.anchors <- FindIntegrationAnchors(rds.list, anchor.features = features)

rna.integrated <- IntegrateData(rds.anchors)
DefaultAssay(rna.integrated) <- "integrated"

rna.integrated <- rna.integrated %>% ScaleData() %>%  
  SCTransform(., assay = "RNA", vars.to.regress = "percent.mt") %>%
  RunPCA(.,assay = "SCT", npcs = NUM_PCA_DIM, verbose = FALSE) %>%
  RunUMAP(., reduction = "pca", reduction.name = "rna.umap", 
          reduction.key = 'rnaUMAP_', dims = 1:NUM_PCA_DIM)
pdf(paste0(savePath, "_UMAP_Integrated_RNA.pdf"))
print(DimPlot(rna.integrated, reduction = "rna.umap", pt.size = 0.4, group.by = "orig.ident"))
dev.off()
cat(sprintf("\n ==>> WTA integration done \n"))

############ Integrate Abseq
rds.list <- lapply(X = rds.list, FUN = function(x) {
  DefaultAssay(x) <- "AB"
  x <- NormalizeData(x, normalization.method = 'CLR', margin = 2)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

features <- SelectIntegrationFeatures(rds.list) # need to run again for ABseq?
rds.list <- lapply(X = rds.list, FUN = function(x) {
  DefaultAssay(x) <- "AB"
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
  return(x)
})
cat(sprintf("\n ==>> After Abseq lapply"))

ab.anchors <- FindIntegrationAnchors(rds.list, assay = rep("AB", length(rds.list)), reduction = "rpca", dims = 1:NUM_ADT_PCA_DIM)
ab.integrated <- IntegrateData(ab.anchors, dims = 1:NUM_ADT_PCA_DIM)
cat(sprintf("\n ==>> After Abseq IntegrateData"))
DefaultAssay(ab.integrated) <- "integrated"

ab.integrated <- ab.integrated %>% ScaleData() %>% RunPCA(., npcs = NUM_ADT_PCA_DIM)
ab.integrated <- RunUMAP(ab.integrated, reduction = "pca", reduction.name = "adt.umap", reduction.key = 'adtUMAP_', dims = 1:NUM_ADT_PCA_DIM)
cat(sprintf("\n ==>> After Abseq UMAP"))

pdf(paste0(savePath, "_UMAP_Integrated_ABSEQ.pdf"))
print(DimPlot(ab.integrated, reduction = "adt.umap", pt.size = 0.4, group.by = "orig.ident"))
dev.off()

## use reciprocal PCA for ADT data integration
sc.obj.anchors <- FindIntegrationAnchors(rds.list, assay = c("AB", "AB"), reduction = "rpca", dims = 1:NUM_ADT_PCA_DIM)
cat(sprintf("\n ==>> After FindIntegrationAnchors"))

# this command creates an 'integrated' data assay
ab.integrated <- IntegrateData(sc.obj.anchors, dims = 1:NUM_ADT_PCA_DIM)
cat(sprintf("\n ==>> After IntegrateData"))

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(ab.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
ab.integrated <- ScaleData(ab.integrated, verbose = FALSE)
cat(sprintf("\n ==>> After ScaleData"))

ab.integrated <- RunPCA(ab.integrated, npcs = NUM_ADT_PCA_DIM, reduction.name = 'apca', verbose = FALSE)
cat(sprintf("\n ==>> After RunPCA"))

ab.integrated <- RunUMAP(ab.integrated, reduction = "apca", reduction.name = "adt.umap", reduction.key = 'adtUMAP_', dims = 1:NUM_ADT_PCA_DIM)
cat(sprintf("\n ==>> After RunUMAP"))

pdf(paste0(project, "_UMAP_Integrated_ADT.pdf"))
print(DimPlot(object = ab.integrated, reduction = "adt.umap", pt.size = 0.4, group.by = "orig.ident"))
dev.off()
DefaultAssau
############ Both for WNN Clustering
rna.integrated[["AB"]] <- ab.integrated
rna.integrated[["apca"]] <- ab.integrated[["apca"]]
rna.integrated[["adt.umap"]] <- ab.integrated[["adt.umap"]]
saveRDS(rna.integrated, file = paste0(savePath, "ADpre_post_C8_44_integrated.RDS"))
cat(sprintf("\n ==>> SaveRDS done \n"))
