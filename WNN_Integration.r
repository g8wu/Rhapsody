#!/share/apps/R/3.4.3/bin/Rscript

library(dplyr)
library(Seurat)
library(ggplot2)
library(EnhancedVolcano)

project <- "Combine_C1_C7"
sink(paste0(project, "_log.txt"))

# number of dimensions from PCA which will be used
# used to plot PCA results, JackStraw plots, and also for clustering
NUM_PCA_DIM <- 30 # 16
NUM_ADT_PCA_DIM <- 18

## Output Directory
OutDir <- getwd()

## read input Seurat objects
sc.obj.list <- list()
obj1 <- readRDS( paste0(OutDir,'/Projects/C6-exact-cutadapt-noIntron_Seurat.rds'))
obj1 <- readRDS( paste0(OutDir,'/Projects/105/105-preAD-exact-cutadapt.rds'))
obj2 <- readRDS( paste0(OutDir,'/Projects/EoE/305/305-preEoE-exact-cutadapt.rds'))

##==============
## first integrate the RNA assay
##==============
sc.obj.list <- list()
sc.obj.list[[1]] <- obj1
sc.obj.list[[2]] <- obj2

# normalize and identify variable features for each dataset independently
sc.obj.list <- lapply(X = sc.obj.list, FUN = function(x) {
		DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = sc.obj.list)

sc.obj.anchors <- FindIntegrationAnchors(object.list = sc.obj.list, anchor.features = features)

# this command creates an 'integrated' data assay
Seurat.combined.RNA <- IntegrateData(anchorset = sc.obj.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Seurat.combined.RNA) <- "integrated"

# Run the standard workflow for visualization and clustering
Seurat.combined.RNA <- ScaleData(Seurat.combined.RNA, verbose = FALSE)
Seurat.combined.RNA <- RunPCA(Seurat.combined.RNA, npcs = NUM_PCA_DIM, verbose = FALSE)
Seurat.combined.RNA <- RunUMAP(Seurat.combined.RNA, reduction = "pca", reduction.name = "rna.umap", reduction.key = 'rnaUMAP_', dims = 1:NUM_PCA_DIM)

pdf(paste0(savePath, "c6_UMAP_Integrated_RNA.pdf"))
print(DimPlot(object = Seurat.combined.RNA, reduction = "rna.umap", pt.size = 0.4, group.by = "orig.ident"))
dev.off()

##==============
## then integrate the ADT assay
##==============

sc.obj.list <- list()
sc.obj.list[[1]] <- ctrl303
sc.obj.list[[2]] <- ctrl306

# normalize and identify variable features for each dataset independently
sc.obj.list <- lapply(X = sc.obj.list, FUN = function(x) {
	DefaultAssay(x) <- "ADT"
  x <- NormalizeData(x, normalization.method = 'CLR', margin = 2)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
cat(sprintf("\n ==>> After lapply"))

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = sc.obj.list)
cat(sprintf("\n ==>> After SelectIntegrationFeatures"))

sc.obj.list <- lapply(X = sc.obj.list, FUN = function(x) {
	DefaultAssay(x) <- "ADT"
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})
cat(sprintf("\n ==>> After lapply - 2"))

## use reciprocal PCA for ADT data integration
sc.obj.anchors <- FindIntegrationAnchors(object.list = sc.obj.list, assay = c("ADT"), reduction = "rpca", dims = 1:NUM_ADT_PCA_DIM)
cat(sprintf("\n ==>> After FindIntegrationAnchors"))

# this command creates an 'integrated' data assay
Seurat.combined.ADT <- IntegrateData(anchorset = sc.obj.anchors, dims = 1:NUM_ADT_PCA_DIM)
cat(sprintf("\n ==>> After IntegrateData"))

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Seurat.combined.ADT) <- "integrated"

# Run the standard workflow for visualization and clustering
Seurat.combined.ADT <- ScaleData(Seurat.combined.ADT, verbose = FALSE)
cat(sprintf("\n ==>> After ScaleData"))

Seurat.combined.ADT <- RunPCA(Seurat.combined.ADT, npcs = NUM_ADT_PCA_DIM, reduction.name = 'apca', verbose = FALSE)
cat(sprintf("\n ==>> After RunPCA"))

Seurat.combined.ADT <- RunUMAP(Seurat.combined.ADT, reduction = "apca", reduction.name = "adt.umap", reduction.key = 'adtUMAP_', dims = 1:NUM_ADT_PCA_DIM)
cat(sprintf("\n ==>> After RunUMAP"))

pdf(paste0(savePath, "ctrl303_306_UMAP_Integrated_ADT.pdf"))
print(DimPlot(object = Seurat.combined.ADT, reduction = "adt.umap", pt.size = 0.4, group.by = "orig.ident"))
dev.off()

##==============
## finally use both of them for clustering + UMAP
##==============

## first add the Seurat.combined.ADT object
## in the ADT assay of the Seurat.combined.RNA
Seurat.combined.RNA[["ADT"]] <- Seurat.combined.ADT[["integrated"]]
Seurat.combined.RNA[["apca"]] <- Seurat.combined.ADT[["apca"]]
Seurat.combined.RNA[["adt.umap"]] <- Seurat.combined.ADT[["adt.umap"]]

## Identify multimodal neighbors. These will be stored in the neighbors slot, 
## and can be accessed using sc.obj[['weighted.nn']]
## The WNN graph can be accessed at sc.obj[["wknn"]], 
## and the SNN graph used for clustering at sc.obj[["wsnn"]]
## Cell-specific modality weights can be accessed at sc.obj$RNA.weight

Seurat.combined.RNA <- FindMultiModalNeighbors(
	Seurat.combined.RNA, 
	reduction.list = list("pca", "apca"), 
	dims.list = list(1:NUM_PCA_DIM, 1:NUM_ADT_PCA_DIM), 
	modality.weight.name = "RNA.weight"
)

##======== UMAP + clustering
Seurat.combined.RNA <- RunUMAP(Seurat.combined.RNA, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

for (res_val in c(0.3, 0.4, 0.5, 0.6, 0.8)) {
	Seurat.combined.RNA <- FindClusters(Seurat.combined.RNA, graph.name = "wsnn", algorithm = 3, resolution = res_val, verbose = FALSE)
	
	##=========== UMAP using WNN
	pdf(paste0(project, "_Clusters_UMAP_WNN_Res", res_val, ".pdf"))
	print(DimPlot(object = Seurat.combined.RNA, reduction = "wnn.umap", pt.size = 0.4))
	print(DimPlot(object = Seurat.combined.RNA, reduction = "wnn.umap", label = TRUE, repel = TRUE, label.size = 6, pt.size = 0.4))
	print(DimPlot(object = Seurat.combined.RNA, reduction = "wnn.umap", pt.size = 0.4, group.by="orig.ident"))
	dev.off()

	## save the object
	SeuratObjFile <- paste0(project, "_SeuratObj_WNN_Res", res_val, ".RDS")
	saveRDS(Seurat.combined.RNA, file= SeuratObjFile)
}

Seurat.combined.RNA <- FindClusters(Seurat.combined.RNA, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)
DimPlot(object = Seurat.combined.RNA, reduction = "wnn.umap", label = TRUE, repel = TRUE, label.size = 6, pt.size = 0.4)
DimPlot(object = Seurat.combined.RNA, reduction = "wnn.umap", pt.size = 0.4, group.by="orig.ident")
SeuratObjFile <- paste0(savePath, "ctrl303_306_SeuratObj_WNN_Res0.8.RDS")
saveRDS(Seurat.combined.RNA, file= SeuratObjFile)

cluster2.markers <- FindMarkers(rds, ident.1 = 2, ident.2 = 4, min.pct = 0.25)
head(cluster2.markers, n = 5)
#clust <- FindAllMarkers(rds, ident.1 = c("1","2"), ident.2 = c("3", "4"), verbose = FALSE)
EnhancedVolcano(cluster2.markers, rownames(cluster2.markers), x = "avg_log2FC", y = "p_val_adj")

volcano <- ggplot(cluster2.markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = ifelse(p_val_adj < 0.05, "Significant", "Not Significant")), alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(x = "Log2(Fold Change)", y = "-log10(P-value)") +
  theme_classic()
