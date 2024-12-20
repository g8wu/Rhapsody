set.seed(99)
library(Seurat)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(tidyverse)
library(EnhancedVolcano)
library(ggplot2)
library(metap)

wkdir <- getwd()

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

anno <- read.csv("coccipilot_anno.csv", header = FALSE)
#rds <- readRDS("NOMID-WNN0.8.rds")
rdsName <- "cocci_pilot_0.5"

####################  Annotate Clusters ####################
annotations <- setNames(anno[, 2], anno[, 1])
rds <- RenameIdents(rds, annotations)

# alphabetize the cell types
Idents(rds) <- factor(Idents(object = rds), levels = sort(levels(rds)))
saveRDS(rds, file = paste0(rdsName, "-annod.RDS"))

write(paste0(Sys.time(), " -> Annotations saved!\n"), file = "LOG.txt", append = TRUE)

# Cell Type Table
write.csv(t(table(Idents(rds))), file=paste0(rdsName, "_cellTypes.csv"))

# PRINT
pdf(paste0(rdsName, "_Annotated.pdf"))
print(DimPlot(rds, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) +
        NoLegend() + plot_annotation(title = paste("Annotated WNN UMAP Louvian: ", rdsName)))
DefaultAssay(rds) <- 'ADT'
annoDotplot <- DotPlot(rds, features = rownames(GetAssayData(rds, assay = 'ADT')),
                       cols = "RdBu", col.min = -1, dot.scale = 8, cluster.idents = TRUE) +
  coord_flip() + labs(title = paste("Annotated AbSeq Dotplot:", rdsName)) + NoLegend()
print(annoDotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
DefaultAssay(rds) <- 'RNA'


DefaultAssay(rds) <- 'SCT'
pdf(paste0(rdsName, "_Annotated.pdf"))
print(DimPlot(rds, reduction = 'rna.umap', label = TRUE, repel = TRUE, label.size = 2.5) +
        NoLegend() + plot_annotation(title = paste("Annotated WNN UMAP Louvian: ", rdsName)))
DefaultAssay(rds) <- 'ADT'
annoDotplot <- DotPlot(rds, features = rownames(GetAssayData(rds, assay = 'ADT')),
                       cols = "RdBu", col.min = -1, dot.scale = 8, cluster.idents = TRUE) +
  coord_flip() + labs(title = paste("Annotated AbSeq Dotplot:", rdsName)) + NoLegend()
print(annoDotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
DefaultAssay(rds) <- 'SCT'

# Unique umi by cell type
metadata <- as.data.frame(FetchData(rds, vars = c("nCount_RNA")))
metadata$annotated_clusters <- as.character(Idents(rds))
readsPerClust <- metadata %>%
  group_by(annotated_clusters) %>%
  summarise(total_reads = sum(nCount_RNA)) %>%
  mutate(annotated_clusters = reorder(annotated_clusters, total_reads))
readsPerClust <- readsPerClust %>%
  mutate(percentage = total_reads / sum(total_reads) * 100)
umiBar <- ggplot(readsPerClust, aes(x = annotated_clusters, y = total_reads)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = paste0(project, ": Total Reads (UMIs) per Annotated Cluster"),
       x = "Cell Type", y = "Total Reads (UMIs)") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = -0.5, size = 4)
print(umiBar)
dev.off()

######### Volcano ##########
n <- length(unique(Idents(rds)))-1
pdf(paste0(project, "_volcanoPlots.pdf"))
for (i in 0:n){
  vector <- setdiff(0:n, i)
  markers <- FindMarkers(rds, ident.1 = i, ident.2 = vector, min.pct = 0.15, verbose = FALSE)
  write.csv(head(markers,), file=paste0(rdsName, "_cellTypes.csv"))
  print(EnhancedVolcano(cluster.markers, rownames(cluster.markers), x = "avg_log2FC", y = "p_val_adj", 
                  title = paste0(project, ": ", i), subtitle = ""))
}  
dev.off()

markers <- FindMarkers(rds, ident.1 = "Eos", min.pct = 0.15, verbose = FALSE)
EnhancedVolcano(markers, rownames(markers), x = "avg_log2FC", y = "p_val_adj", 
                title = "Eos", subtitle = "", )

######### Subcluster ##########
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