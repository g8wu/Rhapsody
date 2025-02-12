# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
## SET WORK DIRECTORY TO SPECIFIC PROJECT FOLDER
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
set.seed(99)
library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(future)
library(sctransform)
library(purrr)
library(leidenAlg)
library(DoubletFinder)
library(plotly)
library(htmlwidgets)

project <- "BALF-Pilot"
fileType <- "mouseOnly-exact-poly10-51x71"
multRate <- 0.03 # Mulitplet rate calculated from all 4 lanes 15.12%
nFeatLower <- 300
nFeatUpper <- 4000
nCountLower <- 150
nCountUpper <-30000
rnaPCs <- 30
adtPCs <- 18
n <- 11     # Let variables reach up to n GB
resRange <- seq(0.7, 1.2, by = 0.1)
clusterAlg <- 4   # FindCluster(alg = 1:  Louvain | fast & effective but not for complex datasets, 
#                                     2:  Refined Louvain | multilevel refined clusters but computationally heavier, 
#                                     3:  Smart Local Moving (SLM) | Louvain w/ more granularity, slower but good for complex datasets
#                                     4:  Leiden | fast & more accurate w/well-connected clusters, requires leidenalg Python package)

##### RENAME ABSEQ
# COCCI
new_names <- c("CD105", "CD115", "CD117", "CD11a", "CD11b", "CD11c", 
               "CD162-Selplg", "CD16/32-FCGR3/2", "CD184", "CD19", "CD31-PECAM1", 
               "CD326-EPCAM", "CD335-NCR1", "CD41-Itga2b", "CD45", "CD48", "CD62L", 
               "CD71-Tfrc", "CXCR2", "Clec7a", "CD150", "F4/80-Adgre1", "Ly6a/E",
               "Ly6G", "NK-1.1-Klrb1b", "SiglecF", "TCRbeta", "TER119/Ly76"
               )

################################################################################
clusterKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")
options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB
wkdir <- getwd()
list <- list()
rds <- readRDS("C:/Users/gio8w/OneDrive - University of California, San Diego Health/Rhapsody/cocci/C60-cocci-mO-exact50k-p10-51_Seurat.rds")
list <- append(list, rds)
# rds_files <- list.files(path = paste0(wkdir,"/"), pattern = paste0("\\-", fileType, ".rds$"), full.names = TRUE)
# list <- map(rds_files, readRDS)

# QC and Rename ADT
for (i in seq_along(list)){
  temp <- GetAssayData(list[[i]], assay = "ADT", layer = "counts")
  rownames(temp) <- new_names
  list[[i]][["ADT"]] <- CreateAssayObject(counts = temp)
  
  list[[i]][["percent.mt"]] <- PercentageFeatureSet(list[[i]], pattern = "^MT-")
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
  list[[i]] <- subset(list[[i]], subset = nFeature_RNA > nFeatLower & nFeature_RNA < nFeatUpper &
                        nCount_RNA > nCountLower & nCount_RNA < nCountUpper)
  # PRINT
  pdf(paste0(wkdir, list[[i]]@meta.data$orig.ident[1], "_QCdVlnPlot.pdf"))
  print(VlnPlot(list[[i]], features = features, pt.size = 0.01, ncol = 3) & labs(x=""))
  dev.off()
}
write(paste0(project, " -> Seurat list built!\n"), file = "LOG", append = TRUE)
gc()

# RNA Normalization -----------------------------------------
for (i in seq_along(list)){
  DefaultAssay(list[[i]]) <- 'RNA'
  #list[[i]] <- SCTransform(list[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = F)
  list[[i]] <- SCTransform(list[[i]], method = "glmGamPoi", verbose = F)   # mouse does not have mitocondrial reads tag
}
gc()
# RNA INTEGRATION -----------------------------------------
# features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
# list <- PrepSCTIntegration(object.list = list, anchor.features = features)
# anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features,
#                                   normalization.method = "SCT")
# combinedRNA <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# gc()
combinedRNA <- list[[1]]
combinedRNA <- RunPCA(combinedRNA, npcs = rnaPCs, assay = "SCT", verbose = FALSE)
combinedRNA <- RunUMAP(combinedRNA, reduction = "pca", reduction.name = "rna.umap",
                       reduction.key = 'rnaUMAP_', dims = 1:rnaPCs)
saveRDS(combinedRNA, file = paste0(project,"-", fileType, "_RNASCT.RDS"))
write(paste0(Sys.time(), " -> RNA SCTransform and integration done and saved!\n"), file = "LOG", append = TRUE)
#PRINT
pdf(paste0(project, "-", fileType, "_RNAUMAP.pdf"))
print(DimPlot(object = combinedRNA, reduction = "rna.umap", pt.size = 0.4, group.by = "orig.ident") + NoLegend())
idents <- unique(combinedRNA$orig.ident)
for (i in seq_along(idents)){
  sample <- subset(combinedRNA, subset = orig.ident == idents[i])
  print(DimPlot(combinedRNA, reduction = "rna.umap", cells.highlight = Cells(sample)) + NoLegend() + 
          labs(title = paste("RNA UMAP", project, ":", idents[i])))
}
dev.off()
gc()

# ADT Normalization -----------------------------------------
for (i in seq_along(list)){
  DefaultAssay(list[[i]]) <- 'ADT'
  list[[i]] <- NormalizeData(list[[i]], normalization.method = 'CLR', margin = 2)
  DefaultAssay(list[[i]]) <- 'RNA'
}
gc()
# ADT Integration -----------------------------------------
# for(i in seq_along(list)){
#   DefaultAssay(list[[i]]) <- 'ADT'
#   list[[i]] <- ScaleData(list[[i]], features = new_names)
#   list[[i]] <- RunPCA(list[[i]], features = new_names, approx = FALSE)
#   DefaultAssay(list[[i]]) <- 'RNA'
# }
# gc()
# assays <- rep(c("ADT"), times = length(list))
# anchors <- FindIntegrationAnchors(object.list = list, assay = assays, 
#                                   reduction = "rpca", dims = 1:adtPCs)
# gc()
# combinedADT <- IntegrateData(anchorset = anchors, dims = 1:adtPCs)
# rm(anchors)
combinedADT <- list[[1]]
combinedADT <- ScaleData(combinedADT, verbose = FALSE)
combinedADT <- FindVariableFeatures(combinedADT)    # FOR SINGLE OBJECT
combinedADT <- RunPCA(combinedADT, npcs = adtPCs, reduction.name = 'apca', verbose = FALSE, approx = FALSE)
combinedADT <- RunUMAP(combinedADT, reduction = "apca", reduction.name = "adt.umap", 
                       reduction.key = 'adtUMAP_', dims = 1:adtPCs)
saveRDS(combinedADT, file = paste0(project, "-", fileType, "_ADTNormalized.RDS"))
write(paste0(Sys.time(), " -> ADT & RNA independent normalization done and saved!\n"), file = "LOG", append = TRUE)
#PRINT
pdf(paste0(project,"-", fileType, "_ADTUMAP.pdf"))
print(DimPlot(object = combinedADT, reduction = "adt.umap", pt.size = 0.4, group.by = "orig.ident"))
dev.off()
write(paste0(Sys.time(), " -> ADT normalization and integration done and saved!\n"), file = "LOG", append = TRUE)

# WNN UMAP WITH RNA & ADT -----------------------------------------
combinedRNA[["ADT"]] <- combinedADT[["ADT"]]
combinedRNA[["apca"]] <- combinedADT[["apca"]]
combinedRNA[["adt.umap"]] <- combinedADT[["adt.umap"]]
saveRDS(combinedRNA, file = paste0(project,"-", fileType, "_integrated.RDS"))
write(paste0(Sys.time(), " -> RNA & ADT INTEGRATION successful, saved before WNN!\n"), file = "LOG", append = TRUE)
rm(combinedADT)
gc()

# DOUBLETFINDER -----------------------------------------
# pK Identification (no ground-truth)
sweep.res.list <- paramSweep(combinedRNA, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
bcmvn <- find.pK(sweep.stats)

ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()

pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

# Homotypic Doublet Proportion Estimate
annotations <- combinedRNA@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(multRate*nrow(combinedRNA@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
combinedRNA <- doubletFinder(combinedRNA, PCs = 1:rnaPCs, pN = 0.25, pK = pK, 
                                nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
View(combinedRNA@meta.data)

#PRINT
pdf(paste0(project,"-", fileType, "_DoubletRNAUMAP.pdf"))
print(DimPlot(combinedRNA, reduction = 'rna.umap', group.by = "DF.classifications_0.25_0.05_863"))
dev.off()
write.csv(table(combinedRNA@meta.data$orig.idents), file=paste0(project, "-", fileType, "_doublets.csv"))
write(paste0(Sys.time(), " -> DoubletFinder done and saved!\n"), file = "LOG", append = TRUE)

combinedRNA <- subset(combinedRNA, subset = DF.classifications_0.25_0.05_4340 == "Singlet")
saveRDS(combinedRNA, file = paste0(project,"-", fileType, "_singlets_integrated.RDS"))

# RNA ONLY WNN -----------------------------------------------------------------
## WNN data: combinedRNA[['weighted.nn']]
## WNN graph: combinedRNA[["wknn"]], 
## SNN graph used for clustering: combinedRNA[["wsnn"]]
## Cell-specific modality weights: combinedRNA$RNA.weight
DefaultAssay(combinedRNA) <- "SCT"
combinedRNA <- FindNeighbors(combinedRNA, dims = 1:rnaPCs, verbose = FALSE)
rds <- combinedRNA
#rds <- RunUMAP(combinedRNA, reduction.name = "rna.umap")
for (res in resRange) {
  rds <- FindClusters(rds, resolution = res, verbose = FALSE, algorithm = clusterAlg)
  saveRDS(rds, file= paste0(project, "-", fileType, "-RNAWNN", res, ".rds"))
  gc()
  #PRINT
  pdf(paste0(project, "-", res, "_RNAWNNUMAP.pdf"))
  print(DimPlot(rds, reduction = 'rna.umap', label = TRUE, repel = TRUE, label.size = 2.5) + 
          NoLegend() + labs(title = paste("RNA WNN UMAP", clusterKey[clusterAlg], res, ":", project)))
  dev.off()
}
write(paste0(Sys.time(), " -> RNA ONLY multi res clustering done! Ready for preview!\n"), file = "LOG", append = TRUE)

# RNA & ADT WNN -----------------------------------------
rds <- FindMultiModalNeighbors(combinedRNA, reduction.list = list("pca", "apca"),
                               dims.list = list(1:rnaPCs, 1:adtPCs), modality.weight.name = "RNA.weight")
rds <- RunUMAP(rds, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
for (res in resRange) {
  rds <- FindClusters(rds, graph.name = "wsnn", algorithm = clusterAlg, resolution = res, verbose = FALSE)
  saveRDS(rds, file= paste0(project, "-ADTRNAWNN", res,  ".rds"))
  gc()
  #PRINT
  pdf(paste0(project, "-", res, "_WNNUMAP.pdf"))
  print(DimPlot(rds, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) +
          NoLegend() + labs(title = paste("ADT&RNA WNN UMAP", clusterKey[clusterAlg], res, ":", project)))
  dev.off()
}
write(paste0(Sys.time(), " -> RNA & ADT multi res clustering done! Ready for preview!\n"), file = "LOG", append = TRUE)

# ______________________________3D umap?
rds <- RunUMAP(rds, reduction = "pca", reduction.name = "3d.umap", 
                                               reduction.key = '3dUMAP_', dims = 1:30, n.components = 3L)
rds <- FindClusters(rds, resolution = 0.5, verbose = FALSE, algorithm = 4)
Embeddings(object = rds, reduction = "3d.umap")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = rds, vars = c("wnnUMAP_1", "wnnUMAP_2", "wnnUMAP_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
fig <- plot_ly(data = plot.data, 
               x = ~wnnUMAP_1, y = ~wnnUMAP_2, z = ~wnnUMAP_3, 
               color = ~seurat_clusters, 
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = .2, width=.5), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text")
saveWidget(fig, paste0(project, "-3DRNAUMAP", res, ".html"))

#PRINT by samples
pdf(paste0(project,"-", fileType, "_RNAUMAP_samples.pdf"))
tags <- unique(rds$Sample_Tag)
tags <- sort(tags)
sample <- c("Multiplet", "C57.1 BALF",  "C57.2 BALF",  "C57.3 BALF",  "DBA.1 BALF", "DBA.2 BALF", "DBA.3 BALF", "Undetermined")
for (tag in tags){
  print(DimPlot(rds, reduction = "rna.umap", cells.highlight = WhichCells(rds, expression = Sample_Tag == tag)) + 
          NoLegend() + labs(title = paste("RNA UMAP", project, sample[1])))
  sample <- sample[-1]
}
dev.off()