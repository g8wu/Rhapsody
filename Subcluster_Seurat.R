#!/share/apps/R/3.4.3/bin/Rscript

library(dplyr)
library(Seurat)
library(ggplot2)
library(MAST)

Input_Seurat_object_filename <- '/mnt/bioadhoc-temp/Groups/vd-ay/Paramita/Single_Cell_RNASeq_DataAnalysis/DATASET/BD_Abseq_Seema_Proj/SourceCode/Compare_AIRR_BDRhap_C1_30k_C7_17k/1_Seurat/Combine_C1_C7_SeuratObj_WNN_Res0.4.RDS'
OUTDIR <- 'out_subclust_7_8_9_13'
project <- 'subclust_7_8_9_13'
ClustLevelVec <- "7:8:9:13"

system(paste("mkdir -p", OUTDIR))

outtextfile <- paste0(OUTDIR, '/', project, '.log')
sink(outtextfile)

cat(sprintf("\n ==>> input parameters - \n Input_Seurat_object_filename : %s \n OUTDIR : %s \n project : %s \n ClustLevelVec : %s ", Input_Seurat_object_filename, OUTDIR, project, ClustLevelVec))

# target cluster levels of cell subset, with respect to the original Seurat object - PARAMITA
target_clust_levels <- as.character(unlist(strsplit(ClustLevelVec, "[,:]")))
cat(sprintf("\n target_clust_levels : %s ", paste(target_clust_levels, collapse=" ")))

#==============
# this R object file stores the Seurat data 
# all the parameters also follow this code
#==============
# number of dimensions from PCA which will be used for clustering
NUM_PCA_DIM <- 30 #16   
NUM_ADT_PCA_DIM <- 18

Target_Res_Val <- 0.4

# read Seurat object
sc.obj <- readRDS(Input_Seurat_object_filename)

# the object has already pre-computed clusters
old_levels <- levels(x = sc.obj)
cat(sprintf("\n ==>> old cluster levels : %s ", paste(as.vector(old_levels), collapse=" : ")))

# now extract only the clusters of "Epithelium"
# basically subset the seurat object
sc.obj.subset <- subset(x = sc.obj, idents = target_clust_levels)

# print new cluster levels (subset)
new_levels <- levels(x = sc.obj.subset)
cat(sprintf("\n ==>> new (subset) cluster levels : %s ", paste(as.vector(new_levels), collapse=" : ")))

#=========================
## Cluster the cells
sc.obj.subset <- FindMultiModalNeighbors(sc.obj.subset, reduction.list = list("pca", "apca"), dims.list = list(1:NUM_PCA_DIM, 1:NUM_ADT_PCA_DIM), modality.weight.name = "RNA.weight")
sc.obj.subset <- RunUMAP(sc.obj.subset, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sc.obj.subset <- FindClusters(sc.obj.subset, graph.name = "wsnn", algorithm = 3, resolution = Target_Res_Val, verbose = FALSE)

pdf(paste0(OUTDIR, '/', project,"_Clusters_UMAP.pdf"))
print(DimPlot(object = sc.obj.subset, reduction = "wnn.umap"))
print(DimPlot(object = sc.obj.subset, reduction = "wnn.umap", label = TRUE, label.size = 4, repel = TRUE))
print(DimPlot(object = sc.obj.subset, reduction = "wnn.umap", pt.size = 0.4, group.by="orig.ident"))
dev.off()

## save current seurat object into a file
## so as to use later
saveRDS(sc.obj.subset, file = paste0(OUTDIR, '/', project, '_single_cell_Seurat_Object.session.rds'))

