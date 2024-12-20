library(Seurat)
library(dplyr)

devtools::install_github("satijalab/seurat-data", force = T)
devtools::install_github("stuart-lab/signac")
devtools::install_github('satijalab/azimuth', ref="develop")

library(Signac)
library(SeuratData)
library(Azimuth)


#c28 <- readRDS("C28-exact_Seurat.rds")
#c28 <- SCTransform(c28) %>% FindVariableFeatures() %>% RunPCA()
#saveRDS(c28, file = "c28_clustered.RDS")
c28 <- readRDS("c28_clustered.RDS")

c28_unds <- c28[,which(c28@meta.data$Sample_Name == "Undetermined")]
c28_unds[["percent.mt"]] <- PercentageFeatureSet(object = c28_unds, pattern = "^MT.")
#c28[["percent.mt"]] <- PercentageFeatureSet(object = c28, pattern = "^MT.")
length(which(c28_unds@meta.data$Sample_Tag == "Undetermined"))


c28_unds <- SCTransform(c28_unds) %>% FindVariableFeatures() %>% RunPCA()
ElbowPlot(c28_unds)

c28_unds <- FindNeighbors(c28_unds, dims = 1:16) %>% RunUMAP(dims = 1:16) %>% FindClusters(resolution = 0.25)
DimPlot(c28_unds)

c28 <- SCTransform(c28) %>% FindVariableFeatures() %>% RunPCA()
ElbowPlot(c28)
c28 <- FindNeighbors(c28, dims = 1:16) %>% RunUMAP(dims = 1:16) %>% FindClusters(resolution = 0.25)

grep("CEACAM8", rownames(c28))

p1 <- VlnPlot(c28_unds, features = "nCount_RNA", group.by = "orig.ident") + ggtitle("nCount_RNA Undetermined Cells")
p2 <- VlnPlot(c28, features = "nCount_RNA", group.by = "orig.ident") + ggtitle("nCount_RNA All Cells")

p1 + p2

p1 <- VlnPlot(c28_unds, features = "percent.mt", group.by = "orig.ident") + ggtitle("percent.mt Undetermined Cells")
p2 <- VlnPlot(c28, features = "percent.mt", group.by = "orig.ident") + ggtitle("percent.mt All Cells")

p1 + p2

p1 <- VlnPlot(c28_unds, features = "CD193-CCR3-AHS0159-pAbO", group.by = "orig.ident") + ggtitle("CCR3 AbSeq Undetermined Cells")
p2 <- VlnPlot(c28, features = "CD193-CCR3-AHS0159-pAbO", group.by = "orig.ident") + ggtitle("CCR3 AbSeq All Cells")

p1 + p2

p1 <- VlnPlot(c28_unds, features = "CD193-CCR3-AHS0159-pAbO", group.by = "orig.ident", y.max = 15) + ggtitle("CCR3 AbSeq Undetermined Cells")
p2 <- VlnPlot(c28, features = "CD193-CCR3-AHS0159-pAbO", group.by = "orig.ident", y.max = 15) + ggtitle("CCR3 AbSeq All Cells")

p1 + p2

p1 <- VlnPlot(c28_unds, features = "SIGLEC8", group.by = "orig.ident", y.max = 15) + ggtitle("SIGLEC8 Undetermined Cells")
p2 <- VlnPlot(c28, features = "SIGLEC8", group.by = "orig.ident", y.max = 15) + ggtitle("SIGLEC8 All Cells")

p1 + p2

p1 <- VlnPlot(c28_unds, features = "CEACAM8", group.by = "orig.ident", y.max = 15) + ggtitle("CEACAM8 Undetermined Cells")
p2 <- VlnPlot(c28, features = "CEACAM8", group.by = "orig.ident", y.max = 15) + ggtitle("CEACAM8 All Cells")

p1 + p2


#saveRDS(c28_unds, file = "/sbgenomics/output-files/c28_undetermined_cells_preprocessed.RDS")

VlnPlot(c28_unds, features = c("CD9", "ITGB9", "PECAM1", "CD36", "ITGA2B", "GP9", "GP1BA", "ITGB3", "PADGEM", "CD63", "CD107a", "Cd40L"), stack = T, flip = T)

c28_unds <- RunAzimuth(c28_unds, ref = "pbmcref")

predictions <- read.delim('/sbgenomics/output-files/azimuth_pred.tsv', row.names = 1)
c28_unds <- AddMetaData(
  object = c28_unds,
  metadata = predictions)

table(c28_unds@meta.data$predicted.celltype.l2)
VlnPlot(c28_unds, features = "predicted.celltype.l2.score", group.by = "predicted.celltype.l2") + ggtitle("Azimuth Confidence Scores Grouped By Cell Type Predictions for Undetermined Cells")




########################################################################################################################################################################
########################################################################################################################################################################
"%notin%" <- Negate("%in%")
st_reads <- read.table(file = "/sbgenomics/project-files/C028/C28_exact/C28-exact_Sample_Tag_ReadsPerCell.csv", header = T, sep = ",")
st_und_cells_reads <- st_reads[which(st_reads$cell %in% colnames(c28_unds)),]
st_det_cells_reads <- st_reads[which(st_reads$cell %notin% colnames(c28_unds)),]
colSums(st_und_cells_reads)
quantile(st_und_cells_reads$SampleTag08_hs_Read_Count)
quantile(st_und_cells_reads$SampleTag09_hs_Read_Count)
quantile(st_und_cells_reads$SampleTag10_hs_Read_Count)

quantile(st_det_cells_reads$SampleTag08_hs_Read_Count)
quantile(st_det_cells_reads$SampleTag09_hs_Read_Count)
quantile(st_det_cells_reads$SampleTag10_hs_Read_Count)

