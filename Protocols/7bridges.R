library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(patchwork)

setwd("C:/Users/Innoscan/Documents/Innoscan_Rhapsody")
cart <- "C036"
savePath <- paste0("RDS_files/",cart, "/")

###################################################################
#rds <- readRDS(file = paste0(savePath, "C36-exact_Seurat.rds"))
rds <- readRDS(paste0(savePath, "clustered.RDS"))
# Integrate datasets
pbmc.big <- merge(pbmc3k, y = c(pbmc4k, pbmc8k), add.cell.ids = c("3K", "4K", "8K"), project = "PBMC15K")
rds <- merge(rds, y = )

# SEPARATE ABSEQ FROM WTA
adt <- rds@assays$RNA@counts[grep("pAbO", rownames(rds)),]
rna <- rds@assays$RNA@counts[-(grep("pAbO", rownames(rds))),]
rds[["ABSEQ"]] <- CreateAssayObject(adt)
abseq <- rownames(rds@assays$ABSEQ)
abseqcat <- c(read.csv("RDS_files/abseqcat.csv", header = F))
rds@assays$ABSEQ@counts@Dimnames[[1]] <- abseqcat$V1

rds[["WTA"]] <- CreateAssayObject(rna)
rds[['RNA']] <- NULL
Assays(rds)
DefaultAssay(rds) = "WTA"

# ============= NORMALIZATION & QC FILTER ============= #
# View QC plots
rds[["percent.mt"]] <- PercentageFeatureSet(rds, pattern = "^MT-")
VlnPlot(rds, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
rds <- subset(rds, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
VlnPlot(rds, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(rds, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(rds, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Normalize Abseq and WTA separately
rds <- NormalizeData(rds, normalization.method = "LogNormalize", scale.factor = 10000)
rds <- NormalizeData(object = rds, normalization.method = "CLR", margin = 2, assay = "ABSEQ") # go across cells, not 

# Find Variable features:
rds <- FindVariableFeatures(rds, selection.method = "vst", nfeatures = 2000)
varPlot <- VariableFeaturePlot(rds)
# Identify the most highly variable genes
top <- head(VariableFeatures(rds), 50)
bot <- tail(VariableFeatures(rds), 50)
varPlot <- LabelPoints(varPlot, points = c(top, bot), repel = TRUE) + title("Overall DEG")
varPlot

write.csv(top, paste0(savePath, paste0(cart,"_top_all_expression.csv")))
write.csv(bot,paste0(savePath, paste0(cart,"bot_all_expression.csv")))

# ==================== CLUSTERING ==================== #
all.genes.rds <- rownames(rds)
rds <- ScaleData(rds)
rds <- RunPCA(rds, features = VariableFeatures(object = rds))
# Examine and visualize PCA results a few different ways
print(rds[["pca"]], dims = 1:20, nfeatures = 5)
VizDimLoadings(rds, dims = 1:2, reduction = "pca")
ElbowPlot(rds)
PCAPlot(rds)

#DimHeatmap(rds, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(rds, dims = 1:30, cells = 500, balanced = TRUE)

############### UMAP CLUSTERING ###############
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
dim <- 20
rds <- JackStraw(rds, num.replicate = 100)
rds <- ScoreJackStraw(rds, dims = 1:dim)
JackStrawPlot(rds, dims = 1:dim)
rds <- FindNeighbors(rds, dims = 1:dim)
rds <- FindClusters(rds, resolution = 0.8)
rds <- RunUMAP(rds, dims = 1:dim)
jpeg(paste0(savePath, cart, "_UMAP.jpeg"), quality = 75, width = 800, height = 800)
UMAPPlot(rds, label = T) + plot_annotation(title = paste("UMAP:", cart)) + theme(legend.position="none")
dev.off()
#rds <- RunTSNE(rds, dims = 1:dim)
#TSNEPlot(rds, label = T) + plot_annotation(title = paste("tSNE: ", saveFileName)) + theme(legend.position="none")

# WNN CLUSTERING
#RunUMAP(rds, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
#DimPlot(object = rds, reduction = "wnn.umap")

# Interactive UMAP
#HoverLocator(plot = umap, information = FetchData(object = rds, vars = c("ident", "PC_1", "nFeature_RNA")))

############### ABSEQ DOTPLOT ###############
# dotplot saved W = 850, H = 750
#rds <- ScaleData(rds)
dotplotTitle = paste0(cart, " AbSeq Dotplot")
jpeg(paste0(savePath,cart, "_Dotplot.jpeg"), quality = 75, width = 1600, height = 800)
DotPlot(rds, features = rownames(GetAssayData(object = rds@assays$ABSEQ, slot = "counts")), 
        cols = "RdBu", col.min = -1, dot.scale = 3, cluster.idents = TRUE) + coord_flip() + plot_annotation(title = paste("AbSeq Dotplot:", cart))
dev.off()
# other cols: Spectral, RdBu

# ABSEQ RIDGEPLOTS
end <- 8
for (i in seq(0, 32, end)) {
  print(i)
  # JPEG Device
  jpeg(paste0(savePath, "ridgeplot_", i, ".jpeg"), quality = 75, width = 2800, height = 400)
  if (i == 24){RidgePlot(rds, assay = "ABSEQ", features = abseq[i:i+7], ncol = end)}
  else{RidgePlot(rds, assay = "ABSEQ", features = abseq[i+1:i+end-1], ncol = end)}
  # Close device
  dev.off()
  print("saved")
}
dev.off()

W = 2800
H = 400
jpeg(paste0(savePath, "ridgeplot_1.jpeg"), quality = 75, width = W, height = H)
RidgePlot(rds, assay = "ABSEQ", features = abseq[0:8], ncol = 8)
dev.off()
jpeg(paste0(savePath, "ridgeplot_2.jpeg"), quality = 75, width = W, height = H)
RidgePlot(rds, assay = "ABSEQ", features = abseq[9:16], ncol = 8)
dev.off()
jpeg(paste0(savePath, "ridgeplot_3.jpeg"), quality = 75, width = W, height = H)
RidgePlot(rds, assay = "ABSEQ", features = abseq[17:24], ncol = 8)
dev.off()
jpeg(paste0(savePath, "ridgeplot_4.jpeg"), quality = 75, width = W, height = H)
RidgePlot(rds, assay = "ABSEQ", features = abseq[25:31], ncol = 8)
dev.off()


############### WTA with ABSEQ Overlay ###############
W = 3000
H = 600
jpeg(paste0(savePath, "featurePlot_1.jpeg"), quality = 75, width = W, height = H)
FeaturePlot(rds, features = abseq[0:5], ncol = 5)
dev.off()
jpeg(paste0(savePath, "featurePlot_2.jpeg"), quality = 75, width = W, height = H)
FeaturePlot(rds, features = abseq[6:10], ncol = 5)
dev.off()
jpeg(paste0(savePath, "featurePlot_3.jpeg"), quality = 75, width = W, height = H)
FeaturePlot(rds, features = abseq[11:15], ncol = 5)
dev.off()
jpeg(paste0(savePath, "featurePlot_4.jpeg"), quality = 75, width = W, height = H)
FeaturePlot(rds, features = abseq[16:20], ncol = 5)
dev.off()
jpeg(paste0(savePath, "featurePlot_5.jpeg"), quality = 75, width = W, height = H)
FeaturePlot(rds, features = abseq[21:25], ncol = 5)
dev.off()
jpeg(paste0(savePath, "featurePlot_6.jpeg"), quality = 75, width = W, height = H)
FeaturePlot(rds, features = abseq[26:30], ncol = 5)
dev.off()
jpeg(paste0(savePath, "featurePlot_7.jpeg"), quality = 75, width = W, height = H)
FeaturePlot(rds, features = abseq[31])
dev.off()

###############  ANNOTATING CLUSTERS ############### 
obj_annotations <- rbind(c("0", "Neutrophil progenitor"),
                         c("1", "Noise"),
                         c("2", "Neutrophil progenitor"),
                         c("3", "Eosinophil"),
                         c("4", "T Cell"),
                         c("5", "Monocyte"),
                         c("6", "T Cell"),
                         c("7", "T Cell"),
                         c("8", "Noise"),
                         c("9", "Neutrophil"),
                         c("10", "Neutrophil"),
                         c("11", "T Cell"),
                         c("12", "Natural Killer"),
                         c("13", "Platelet"),
                         c("14", "B cell"),
                         c("15", "Monocyte"),
                         c("16", "Activated B cell"),
                         c("17", "Basophil"),
                         c("18", "HPSC"))
colnames(obj_annotations) <- c("Cluster", "CellType")
obj_annotations <- data.frame(obj_annotations)

# save the annotations as a csv if you'd like
# write.csv(obj_annotations, file = file.path(path_data, "obj_annotations.csv"))

# prepare the annotation information
annotations <- setNames(obj_annotations$CellType, obj_annotations$Cluster)

# relabel the Seurat clusters
# Idents(obj) <- "seurat_clusters"
rds <- RenameIdents(object = rds, annotations)

# alphabetize the cell types
Idents(rds) <- factor(Idents(object = rds), levels = sort(levels(rds)))

# useful metadata (e.g. if you want to have multiple annotation sets)
rds[["annotated_clusters"]] <- Idents(object = rds)

# save the processed and annotated Seurat object
# saveRDS(obj, file = file.path(path_data, "tenx_pbmc5k_CITEseq_annotated.rds"))

# info about the clusters
obj_annotations %>%
  group_by(CellType) %>%
  transmute(Clusters = paste0(Cluster, collapse = ", ")) %>%
  distinct() %>%
  arrange(CellType)


jpeg(paste0(savePath, cart, "_annotatedUMAP.jpeg"), quality = 75, width = 800, height = 800)
UMAPPlot(rds, label = T) + plot_annotation(title = paste("Annotated UMAP:", cart)) + theme(legend.position="none")
dev.off()
## UMAP Cells/cluster
write.csv(t(table(Idents(rds))), file=paste0(savePath, paste0(cart, "_clusters.csv")))

DotPlot(rds, features = rownames(GetAssayData(object = rds@assays$ABSEQ, slot = "counts")), 
        cols = "RdBu", col.min = -1, dot.scale = 3, cluster.idents = TRUE) + coord_flip() + RotatedAxis() + plot_annotation(title = paste("AbSeq Dotplot:", cart))

# metadata file for MCIA
meta <- data.frame("CellType" = Idents(rds))




############### SAVE RDS ###############
DefaultAssay(rds) = "WTA" 
saveRDS(rds, paste0(savePath, "clustered.RDS"))
