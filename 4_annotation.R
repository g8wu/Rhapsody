install.packages('devtools')
devtools::install_github('immunogenomics/presto')
BiocManager::install('multtest')
set.seed(99)
library(Seurat)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(tidyverse)
if(!require("EnhancedVolcano", quietly = T)) BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
if(!require("metap", quietly = T)) install.packages("metap")
library(metap)

res <- "0.3"
project <- "BRAIN"

# Get command line arguments
# args <- commandArgs(trailingOnly = TRUE)
# 
# # Check if arguments are provided
# if (length(args) == 0) {
#   stop("No files provided. Please input: \n1 - clustered Seurat object for annotation \n2 - annotation .csv file.")
# }
# print(args[1])
# print(args[2])
# # Excel file with Col 1: cluster number, Col 2: Cell type annotation
# anno <- read.csv(paste0(wkdir, "/", args[1]))
# rds <- readRDS(paste0(wkdir, "/", args[2]))
# rdsName <- substr(args[1], 1, nchar(args[1]) - 4)

####################  Annotate Clusters ####################
anno <- read.csv(paste(res, project, "anno.csv", sep  ="-"), header = FALSE)
annotations <- setNames(anno[, 2], anno[, 1])
rds <- RenameIdents(rds, annotations)
rds@meta.data$annotations <- Idents(rds)
table(rds$annotations)

# alphabetize the cell types
Idents(rds) <- factor(Idents(object = rds), levels = sort(levels(rds)))

# Cell Type Table
write.csv(table(Idents(rds)), file=paste(project, "ANNO-cellTypes.csv", sep = "-"))

# PRINT UMAP
pdf(paste(project, res, "ANNO-RNAClust.pdf", sep = "-"))
print(DimPlot(rds, reduction = 'rna.umap', label = TRUE, repel = TRUE, label.size = 2.5) +
        NoLegend() + plot_annotation(title = paste("Annotated RNA UMAP Leiden: ", project)))
dev.off()

# PRINT ANNO Abseq Dotplot
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
colors <- colorRampPalette(brewer.pal(n = 9, name = "Oranges"))
pdf(paste(project, res, "ANNO-AbseqDot.pdf", sep = "-"),width = 15, height = 10)
Clustered_DotPlot(rds, features = abseq, colors_use_exp = colors(20), 
                  exp_color_min = 0, show_ident_legend = F)
dev.off()
DefaultAssay(rds) <- "SCT"

# QC post Anno 
features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_ADT")

pdf(paste(project, res, "ANNO-VlnPlot.pdf", sep = "-"), width = 20, height = 6)
print(VlnPlot(rds, features = features, pt.size = 0, ncol = length(features), group.by = "annotations"))
dev.off()

############# GENES Dot Plot Annod ############# 
rds$orig.clust = paste0(rds$annotations, "_", rds$orig.ident )

# START!!
DefaultAssay(rds) <- "SCT"
genes <- read.csv("CellType.csv", header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)
pdf(paste(project, res, "ANNO-CellTypeList-SCTDot.pdf", sep = "-"))
for (col in names){
  print(col)
  print(DotPlot(rds, features = genes[[col]], col.min = 0, group.by = "orig.clust") + 
          labs(title = paste("SCT Dotplot:", col)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}
dev.off()

# START!!
DefaultAssay(rds) <- "SCT"
genes <- read.csv("InterestGenes.csv", header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)
pdf(paste(project, res, "ANNO-InterestList-SCTDot.pdf", sep = "-"), height = 11, width = 13)
for (col in names){
  print(col)
  print(DotPlot(rds, features = genes[[col]], col.min = 0, group.by = "orig.clust") + 
          labs(title = paste("SCT Dotplot:", col)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}
dev.off()

######### Cluster DEG & Volcano POST ANNO ##########
project = "SPLEEN"
table(rds$annotations)
rds <- PrepSCTFindMarkers(rds)
degs <- FindAllMarkers(rds, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(degs, file=paste(project, res, "ANNO-DEGs.csv", sep = "-"))
saveRDS(rds, file = paste(project, res, "ANNO-DEGs.RDS", sep = "-"))

pdf(paste(project, res, "ANNO-Volcano.pdf", sep = "-"))
for (i in rownames(table(Idents(rds)))){
  clustDegs <- subset(degs, cluster == i)
  print(EnhancedVolcano(clustDegs, rownames(clustDegs), x = "avg_log2FC", y = "p_val_adj",
                        title = paste(project, "cluster", i, sep = " "), subtitle = ""))
}
dev.off()

# Test with no ad p val == 0
pdf(paste(project, res, "ANNO-Volcano-NoAdjP0.pdf", sep = "-"))
for (i in rownames(table(Idents(rds)))){
  clustDegs <- subset(degs, cluster == i)
  clustDegs <- clustDegs[clustDegs$p_val_adj != 0,]   # for taking out adj p val == 0
  print(EnhancedVolcano(clustDegs, rownames(clustDegs), x = "avg_log2FC", y = "p_val_adj",
                        title = paste(project, "cluster", i, sep = " "), subtitle = ""))
}
dev.off()

project = "BRAIN"
rds <- readRDS("~/nomid/BRAIN/BRAIN-nomid-exact-p10-51x71-0.3-RNAClust.rds")
anno <- read.csv(paste(res, project, "anno.csv", sep  ="-"), header = FALSE)
annotations <- setNames(anno[, 2], anno[, 1])
rds <- RenameIdents(rds, annotations)
rds@meta.data$annotations <- Idents(rds)
table(rds$annotations)

# alphabetize the cell types
Idents(rds) <- factor(Idents(object = rds), levels = sort(levels(rds)))
table(rds$annotations)
rds <- PrepSCTFindMarkers(rds)
degs <- FindAllMarkers(rds, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(degs, file=paste(project, res, "ANNO-DEGs.csv", sep = "-"))
saveRDS(rds, file = paste(project, res, "ANNO-DEGs.RDS", sep = "-"))

pdf(paste(project, res, "ANNO-Volcano.pdf", sep = "-"))
for (i in 1:length(unique(rds$annotations))){
  clustDegs <- subset(degs, annotations == i)
  print(EnhancedVolcano(clustDegs, rownames(clustDegs), x = "avg_log2FC", y = "p_val_adj",
                        title = paste(project, "cluster", i, sep = " "), subtitle = ""))
}
dev.off()

# Test with no ad p val == 0
pdf(paste(project, res, "ANNO-Volcano-NoAdjP0.pdf", sep = "-"))
for (i in 1:length(unique(rds$annotations))){
  clustDegs <- subset(degs, annotations == i)
  clustDegs <- clustDegs[clustDegs$p_val_adj != 0,]   # for taking out adj p val == 0
  print(EnhancedVolcano(clustDegs, rownames(clustDegs), x = "avg_log2FC", y = "p_val_adj",
                        title = paste(project, "cluster", i, sep = " "), subtitle = ""))
}
dev.off()

######### WT VS ND DEG & Volcano POST ANNO ##########
project = "SPLEEN"
rds <- readRDS("/SPLEEN/SPLEEN-nomid-exact-p10-51x71-0.3-RNAClust.rds")
table(rds$annotations)
Idents(rds) <- rds$orig.idents
rds <- PrepSCTFindMarkers(rds)
pdf(paste(project, res, "ANNO-WTvNMDVolcano.pdf", sep = "-"))
for (i in 1:length(unique(rds$annotations))){
  clust <- subset(rds, annotations == i)
  clustDegs <- FindMarkers(clustDegs, ident.1 =  paste0("NOMID-", project), ident.2 = paste0("WT-", project), 
                           only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  clustDegs <- arrange(clustDegs, desc(avg_log2FC))
  write.csv(clustDegs, file=paste(project, res, "ANNO-DEGs-WTvNMD.csv", sep = "-"))
  print(EnhancedVolcano(clustDegs, rownames(clustDegs), x = "avg_log2FC", y = "p_val_adj",
                        title = paste(project, "cluster", i, sep = " "), subtitle = ""))
}
dev.off()

project = "BRAIN"
rds <- readRDS("/BRAIN/BRAIN-nomid-exact-p10-51x71-0.3-RNAClust.rds")
table(rds$annotations)
Idents(rds) <- rds$orig.idents
rds <- PrepSCTFindMarkers(rds)
pdf(paste(project, res, "ANNO-WTvNMDVolcano.pdf", sep = "-"))
for (i in 1:length(unique(rds$annotations))){
  clust <- subset(rds, annotations == i)
  clustDegs <- FindMarkers(clustDegs, ident.1 =  paste0("NOMID-", project), ident.2 = paste0("WT-", project), 
                           only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
  clustDegs <- arrange(clustDegs, desc(avg_log2FC))
  write.csv(clustDegs, file=paste(project, res, "ANNO-DEGs-WTvNMD.csv", sep = "-"))
  print(EnhancedVolcano(clustDegs, rownames(clustDegs), x = "avg_log2FC", y = "p_val_adj",
                        title = paste(project, "cluster", i, sep = " "), subtitle = ""))
}
dev.off()

############## Subcluster Annod ############## 
clusterKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")

cluster <- "13"
res <- 0.3
clusterAlg <- 2
ident1 <- 0
ident2 <- 1
select <- WhichCells(rds, idents = cluster)
sub <- subset(rds, cells = select)
sub <- RunPCA(sub, assay = "SCT", verbose = FALSE)
DefaultAssay(sub) <-"SCT"
DimPlot(sub, reduction = "pca")
sub <- RunUMAP(sub, reduction = "pca", reduction.name = "rna.umap",
               reduction.key = 'rnaUMAP_', dims = 1:50)
DimPlot(sub, reduction = "rna.umap")
sub <- FindClusters(sub, resolution = res, verbose = FALSE, algorithm = clusterAlg)
DimPlot(sub, reduction = 'rna.umap', label = TRUE, repel = TRUE, label.size = 5) + 
  NoLegend() + labs(title = paste(cluster, "RNA WNN UMAP", res, clusterKey[clusterAlg]))
markers <- FindMarkers(sub, ident.1 = ident1, ident.2 = ident2, min.pct = 0.15, verbose = FALSE)
EnhancedVolcano(markers, rownames(markers), x = "avg_log2FC", y = "p_val_adj", 
                title = paste("Cluster ", cluster,"Subcluster", ident1, "v", ident2), subtitle = "", ) + NoLegend()

tags <- unique(sub$Sample_Tag)
tags <- sort(tags)
sample <- c("Multiplet", "C57.1 BALF",  "C57.2 BALF",  "C57.3 BALF",  "DBA.1 BALF", "DBA.2 BALF", "DBA.3 BALF", "Undetermined")

for (tag in tags){
  png(paste0("Cluster ", cluster, " subcluster ", sample[1], ".png"), width = 800, height = 800)
  print(DimPlot(sub, sizes.highlight = 3, reduction = "rna.umap", cells.highlight = WhichCells(sub, expression = Sample_Tag == tag)) + 
          NoLegend() + labs(title = paste("RNA UMAP", sample[1])))
  sample <- sample[-1]
  dev.off()
}


############## CELLCHAT ######################

library(CellChat)
library(patchwork)
library(circlize)

# Create CellChat object
DefaultAssay(rds) <- "SCT"
data.input <- GetAssayData(rds, assay = "SCT", slot = "data")
labels <- Idents(rds)
meta <- data.frame(group = labels, row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# Set CellChat database
CellChatDB <- CellChatDB.human  # Use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

# Subset data
cellchat <- subsetData(cellchat)

# Identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Optional: project gene expression data onto protein-protein interaction (PPI)
cellchat<- projectData(cellchat, PPI.human)

# Compute communication probability
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)

# Filter out communications if low cell numbers
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Save the CellChat object
saveRDS(cellchat, file = "EoE-312_12_14prepostexact-noIntron_annotatedWNN0.8cellchat.rds")
#cellchat <- readRDS("C:/Users/gio8w/OneDrive - University of California, San Diego Health/Rhapsody/EoE/EoE-312_12_14prepostexact-noIntron_annotatedWNN0.8cellchat.rds")

cellchat <- aggregateNet(cellchat)
cellchat@net$count
cellchat@net$weight
available_signaling <- unique(cellchat@netP$pathways)
write.csv(available_signaling, file = "cellChat signals.csv")

# Visualize interaction networks
groupSize <- as.numeric(table(cellchat@idents))
idents <- unique(cellchat@idents)
png("CellChat_EoEprepost302_312_214.png", width = 800, height = 800)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, vertex.size = 5, 
                 vertex.label.cex = 2,edge.width = 0.5,title.name = "EoE pre/post 301, 312, 314")
dev.off()

png("CellChat_typeEoEprepost302_312_214.png", width = 2000, height = 1000)
par(mfrow = c(2, 5), xpd=T)
for (i in idents){
  netVisual_circle(cellchat@net$count, idents.use = c(i),
                   vertex.weight = groupSize, vertex.size = 5, 
                   vertex.label.cex = 2,edge.width = 0.5, title.name = i)
}
mtext("Cell-Cell Communication Network: Uninfected WT BALF", side = 3, outer = T, 
      line = -2)
dev.off()
netVisual_bubble(cellchat, title = "Bubble Plot Signaling")
# Visualize interaction heatmap
netVisual_heatmap(cellchat)



signals <- c("PECAM1", "IL16", "CLEC", "ITGB2", "TGFb")
for (signal in signals){
  png(paste0("CellChat path ", signal, ".png"))
  # Visualize signaling pathways
  netVisual_aggregate(cellchat, signaling = signal, layout = "circle")
  title(main = signal)
  dev.off()
}

# Unique umi by cell type
# metadata <- as.data.frame(FetchData(rds, vars = c("nCount_RNA")))
# metadata$annotated_clusters <- as.character(Idents(rds))
# readsPerClust <- metadata %>%
#   group_by(annotated_clusters) %>%
#   summarise(total_reads = sum(nCount_RNA)) %>%
#   mutate(annotated_clusters = reorder(annotated_clusters, total_reads))
# readsPerClust <- readsPerClust %>%
#   mutate(percentage = total_reads / sum(total_reads) * 100)
# umiBar <- ggplot(readsPerClust, aes(x = annotated_clusters, y = total_reads)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   labs(title = paste0(project, ": Total Reads (UMIs) per Annotated Cluster"),
#        x = "Cell Type", y = "Total Reads (UMIs)") +
#   geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = -0.5, size = 4)
# pdf(paste(project, res, "ANNO-UMIBarplot.pdf", sep = "-"))
# print(umiBar)
# dev.off()