# DoubletFinder
if(!require("remotes", quietly = T)) install.packages("remotes")
if (!require("DoubletFinder",quietly = T)) remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

library(DoubletFinder)
library(Seurat)
library(ggplot2)
library(purrr)
library(dplyr)
set.seed(99)

project <- "dress"
fileType <- "exact-p10-51x71"
multRate <- 0.05 
rnaPCs <- 30
adtPCs <- 20
n <- 40     # Let variables reach up to n GB
options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB
resRange <- seq(0.3, 0.4, by = 0.1)
# FindCluster(alg = 1:  Louvain | fast & effective but not for complex datasets,
#                   2:  Refined Louvain | multilevel refined clusters but computationally heavier,
#                   3:  Smart Local Moving (SLM) | Louvain w/ more granularity, slower but good for complex datasets
#                   4:  Leiden | fast & more accurate w/well-connected clusters, requires leidenalg Python package)
clusterAlg <- 4
clusterKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")
wkdir <- getwd()

rds <- readRDS("BRAIN-ANNO.RDS")
print("Seurat object read in---------------------------")
rds
print("------------------------------------------------")

# DOUBLETFINDER -----------------------------------------
# pK Identification (no ground-truth)
sweep.res.list <- paramSweep(rds, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()

pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

# Homotypic Doublet Proportion Estimate
annotations <- rds@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(multRate*nrow(rds@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder
rds <- doubletFinder_v3(rds, PCs = 1:rnaPCs, pN = 0.25, pK = pK,
                                nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
#PRINT
pdf(paste0(project,"-", fileType, "_DoubletRNAUMAP.pdf"))
print(DimPlot(rds, reduction = 'umap'))
print(DimPlot(rds, reduction = 'umap', group.by = "orig.ident"))
dev.off()

write.csv(table(rds@meta.data$orig.idents), file=paste0(rds@project.name, "_doublets.csv"))
saveRDS(rds, file = paste0(rds@project.name, "_integrated_Doublets.RDS"))
print(paste0(Sys.time(), " -> DoubletFinder done and saved!\n"))