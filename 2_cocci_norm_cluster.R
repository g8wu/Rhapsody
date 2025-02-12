# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
# ## SET WORK DIRECTORY TO SPECIFIC PROJECT FOLDER
# ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
# system("pip install leidenalg")
# py_install("pandas")
# set.seed(99)
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if(!require("remotes", quietly = T)) install.packages("remotes")
# library(Seurat)
# if (!require("harmony", quietly = TRUE)) BiocManager::install("harmony")
# library(harmony)
# library(tidyverse)
# library(sctransform)
# BiocManager::install("glmGamPoi")
# if (!require("glmGamPoi", quietly = T)) BiocManager::install("glmGamPoi")
# library(glmGamPoi)
# library(ggplot2)
# if(!require("future", quietly = T)) install.packages("future")
# library(future)
# library(sctransform)
# library(purrr)
# if(!require("leidenAlg", quietly= T)) install.packages("leidenAlg")
# library(leidenAlg)
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
# library(DoubletFinder)
# if(!require("leidenAlg", quietly= T)) install.packages("reticulate")
# library(reticulate)


project <- "cocci-"
fileType <- "exact-p10-51x71"
multRate <- 0.05 # Mulitplet rate calculated from all 4 lanes 15.12%
rnaPCs <- 20
adtPCs <- 20
n <- 40     # Let variables reach up to n GB
options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB
resRange <- seq(0.3, 0.6, by = 0.1)
# FindCluster(alg = 1:  Louvain | fast & effective but not for complex datasets,
#                   2:  Refined Louvain | multilevel refined clusters but computationally heavier,
#                   3:  Smart Local Moving (SLM) | Louvain w/ more granularity, slower but good for complex datasets
#                   4:  Leiden | fast & more accurate w/well-connected clusters, requires leidenalg Python package)
clusterAlg <- 4
clusterKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")
wkdir <- getwd()


##### RENAME ABSEQ
# COCCI
#rownames(rds@assays$ADT@data)
new_names <- c("CD105", "CD115", "CD117", "CD11a", "CD11b", "CD11c", "CD162", "CD16/32", "CD184", "CD19", "CD31",
               "CD326", "CD335", "CD41", "CD45", "CD48", "CD62L", "CD71", "CXCR2", "Clec7a", "CD150", "F4/80", "Ly6A/E",
               "Ly6G", "NK1.1", "SiglecF", "TCRb", "TER119")
rds_files <- list.files(path = paste0(wkdir, "/samples"), pattern = paste0("\\-mouseCocci-LUNG-exact50k.rds$"), full.names = TRUE)
print(rds_files)
list <- list()
list <- map(rds_files, readRDS)

# QC and Rename ADT
pdf(paste0(project, "-", fileType, "_VlnPlot.pdf"))
#features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cocci")
features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
for (i in seq_along(list)){
  temp <- list[[i]]@assays$ADT@counts
  rownames(temp) <- new_names
  list[[i]][["ADT"]] <- CreateAssayObject(counts = temp)
  rm(temp)
  list[[i]][["percent.mt"]] <- PercentageFeatureSet(list[[i]], pattern = "mt-")
  #list[[i]][["percent.cocci"]] <- PercentageFeatureSet(list[[i]], pattern = "CIMG")
  print(VlnPlot(list[[i]], features = features, pt.size = 0.01, ncol = length(features), same.y.lims = T) & 
    theme(axis.text.x = element_blank()) & labs(list[[i]]$orig.ident))
}
dev.off()

# Set QC thresholds
nFeatLower <- 300
nFeatUpper <- 4000
nCountLower <- 150
nCountUpper <-15000
mtPct <- 20
features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
pdf(paste0(project, "-", fileType, "_QCdVlnPlot.pdf"))
for (i in seq_along(list)){
  list[[i]] <- subset(list[[i]], subset = nFeature_RNA > nFeatLower & nFeature_RNA < nFeatUpper &
                        nCount_RNA > nCountLower & nCount_RNA < nCountUpper)
  print(VlnPlot(list[[i]], features = features, pt.size = 0.01, ncol = 4) & 
          theme(axis.text.x = element_blank()) & labs(x=list[[i]]$orig.ident))
}
dev.off()
write(paste0(project, " -> Seurat list QCd and built!\n"), file = "LOG.txt", append = TRUE)
gc()

# RNA Normalization -----------------------------------------
for (i in seq_along(list)){
  DefaultAssay(list[[i]]) <- 'RNA'
  list[[i]] <- SCTransform(list[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = F)
}
gc()

# RNA INTEGRATION -----------------------------------------
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 2000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features,
                                  normalization.method = "SCT")
combinedRNA <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
gc()
combinedRNA <- RunPCA(combinedRNA, npcs = rnaPCs, verbose = FALSE)
combinedRNA <- RunUMAP(combinedRNA, reduction = "pca", reduction.name = "rna.umap",
                       reduction.key = 'rnaUMAP_', dims = 1:rnaPCs)
saveRDS(combinedRNA, file = paste0(project,"-", fileType, "_RNASCT.RDS"))
write(paste0(Sys.time(), " -> RNA SCTransform and integration done and saved!\n"), file = "LOG.txt", append = TRUE)
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
for(i in seq_along(list)){
  DefaultAssay(list[[i]]) <- 'ADT'
  list[[i]] <- ScaleData(list[[i]], features = new_names)
  list[[i]] <- RunPCA(list[[i]], features = new_names, approx = FALSE)
  DefaultAssay(list[[i]]) <- 'RNA'
}
gc()
assays <- rep(c("ADT"), times = length(list))
anchors <- FindIntegrationAnchors(object.list = list, assay = assays,
                                  reduction = "rpca", dims = 1:adtPCs)
gc()
combinedADT <- IntegrateData(anchorset = anchors, dims = 1:adtPCs)
rm(anchors)
combinedADT <- ScaleData(combinedADT, verbose = FALSE)
combinedADT <- RunPCA(combinedADT, npcs = adtPCs, reduction.name = 'apca', verbose = FALSE, approx = FALSE)
combinedADT <- RunUMAP(combinedADT, reduction = "apca", reduction.name = "adt.umap",
                       reduction.key = 'adtUMAP_', dims = 1:adtPCs)
saveRDS(combinedADT, file = paste0(project, "-", fileType, "_ADTNormalized.RDS"))
write(paste0(Sys.time(), " -> ADT & RNA independent normalization done and saved!\n"), file = "LOG.txt", append = TRUE)
#PRINT
pdf(paste0(project,"-", fileType, "_ADTUMAP.pdf"))
print(DimPlot(object = combinedADT, reduction = "adt.umap", pt.size = 0.4, group.by = "orig.ident"))
dev.off()
write(paste0(Sys.time(), " -> ADT normalization and integration done and saved!\n"), file = "LOG.txt", append = TRUE)

# INTEGRATION -----------------------------------------
combinedRNA[["ADT"]] <- combinedADT[["ADT"]]
combinedRNA[["apca"]] <- combinedADT[["apca"]]
combinedRNA[["adt.umap"]] <- combinedADT[["adt.umap"]]

# Batch effect correction
DefaultAssay(combinedRNA) <- "RNA"

# run on RNA and store
combinedRNA <- RunHarmony(combinedRNA, group.by.vars = "orig.ident", assay.use = "RNA")
harmony <- Embeddings(combinedRNA, "harmony")
combinedRNA[["harmony_RNA"]] <- CreateDimReducObject(harmony, key = "harmonyRNA_", assay = "RNA")

# run on ADT and store
combinedRNA <- RunHarmony(combinedRNA, group.by.vars = "orig.ident", assay.use = "ADT")
harmony <- Embeddings(combinedRNA, "harmony")
combinedRNA[["harmony_ADT"]] <- CreateDimReducObject(harmony, key = "harmonyADT_", assay = "ADT")

# PRINT results
pdf(paste0(project, "-harmony.pdf"))
DimPlot(combinedRNA, reduction = "harmony_RNA", group.by = "orig.ident")
DimPlot(combinedRNA, reduction = "harmony_ADT", group.by = "orig.ident")
dev.off()

# SAVE!
saveRDS(combinedRNA, file = paste0(project,"-", fileType, "_integrated.RDS"))
write(paste0(Sys.time(), " -> RNA & ADT INTEGRATION successful, saved before WNN!\n"), file = "LOG.txt", append = TRUE)
rm(combinedADT)
rm(harmony)
gc()


# DOUBLETFINDER -----------------------------------------
# pK Identification (no ground-truth)
# sweep.res.list <- paramSweep_v3(combinedRNA, PCs = 1:20, sct = FALSE)
# sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
# bcmvn <- find.pK(sweep.stats)
#
# ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
#
# pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
#   filter(BCmetric == max(BCmetric)) %>%
#   select(pK)
# pK <- as.numeric(as.character(pK[[1]]))
#
# # Homotypic Doublet Proportion Estimate
# annotations <- combinedRNA@meta.data$seurat_clusters
# homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
# nExp_poi <- round(multRate*nrow(combinedRNA@meta.data))
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#
# # run doubletFinder
# combinedRNA <- doubletFinder_v3(combinedRNA, PCs = 1:rnaPCs, pN = 0.25, pK = pK,
#                                 nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
#
# #PRINT
# pdf(paste0(project,"-", fileType, "_DoubletRNAUMAP.pdf"))
# print(DimPlot(combinedRNA, reduction = 'umap'))
# print(DimPlot(combinedRNA, reduction = 'umap', group.by = "orig.ident"))
# dev.off()
# write.csv(table(combinedRNA@meta.data$orig.idents), file=paste0(project, "-", fileType, "_doublets.csv"))
# saveRDS(combinedRNA, file = paste0(project,"-", fileType, "_integrated_Doublets.RDS"))
# write(paste0(Sys.time(), " -> DoubletFinder done and saved!\n"), file = "LOG.txt", append = TRUE)

# RNA ONLY WNN -----------------------------------------------------------------
## WNN data: combinedRNA[['weighted.nn']]
## WNN graph: combinedRNA[["wknn"]],
## SNN graph used for clustering: combinedRNA[["wsnn"]]
## Cell-specific modality weights: combinedRNA$RNA.weight
DefaultAssay(combinedRNA) <- "RNA"
rds <- FindVariableFeatures(combinedRNA, assay = "RNA", selection.method = "vst", nfeatures = 2000)
rds <- FindNeighbors(rds, dims = 1:rnaPCs, graph.name = "RNA_nn")
rds <- RunUMAP(rds, dims = 1:rnaPCs)
for (res in resRange) {
  rds <- FindClusters(rds, algorithm = clusterAlg, resolution = res, graph.name = "RNA_nn", verbose = FALSE)
  saveRDS(rds, file= paste0(project, "-", fileType, "-RNAClust", res, ".rds"))
  gc()
  #PRINT
  pdf(paste0(project, "-", res, "_RNAWNNUMAP.pdf"))
  print(DimPlot(rds, reduction = 'rna.umap', label = TRUE, repel = TRUE, label.size = 2.5) +
          NoLegend() + labs(title = paste("RNA WNN UMAP", clusterKey[clusterAlg], res, ":", project)))
  dev.off()
}
write(paste0(Sys.time(), " -> RNA ONLY multi res clustering done! Ready for preview!\n"), file = "LOG.txt", append = TRUE)

# RNA & ADT WNN -----------------------------------------
rds <- FindVariableFeatures(rds, assay = "RNA", selection.method = "vst", nfeatures = 500)
rds <- FindMultiModalNeighbors(rds, reduction.list = list("pca", "apca"),
                               dims.list = list(1:rnaPCs, 1:adtPCs), modality.weight.name = "integrated.weight")
rds <- RunUMAP(rds, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
for (res in resRange) {
  rds <- FindClusters(rds, graph.name = "wsnn", algorithm = clusterAlg, resolution = res, verbose = FALSE)
  saveRDS(rds, file= paste0(project, "-WNN", res, ".rds"))
  gc()
  #PRINT
  pdf(paste0(project, "-", res, "_WNNUMAP.pdf"))
  print(DimPlot(rds, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) +
          NoLegend() + labs(title = paste("ADT&RNA WNN UMAP", clusterKey[clusterAlg], res, ":", project)))
  dev.off()
}
write(paste0(Sys.time(), " -> RNA & ADT multi res clustering done! Ready for preview!\n"), file = "LOG.txt", append = TRUE)
