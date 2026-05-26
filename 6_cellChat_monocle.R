# CellChat ######################
# Tutorial: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#part-i-data-input-processing-and-initialization-of-cellchat-object
# Manual: https://www.rdocumentation.org/packages/CellChat/versions/1.0.0
if(!require("ComplexHeatmap", quietly = T)) BiocManager::install("ComplexHeatmap")
if (!require("CellChat", quietly = TRUE)) devtools::install_github("jinworks/CellChat")
library(CellChat)
library(patchwork)
library(circlize)

# Reorder CellChat idents
cellchat@idents <- factor(cellchat@idents, levels = sort(levels(cellchat@idents)))

## Create CellChat object ####
group <- ""
Idents(rds) <- group
DefaultAssay(rds) <- "SCT"
data.input <- GetAssayData(rds, assay = "SCT", slot = "data")
labels <- Idents(rds)
meta <- data.frame(group = labels, row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# Set CellChat database
# CellChatDB <- CellChatDB.mouse  # mouse data
CellChatDB <- CellChatDB.human  # human data
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

# Subset data
cellchat <- subsetData(cellchat)

# Identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Optional: project gene expression data onto protein-protein interaction (PPI)
# cellchat<- projectData(cellchat, PPI.mouse)
cellchat<- projectData(cellchat, PPI.human)

# Compute communication probability
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)

# Filter out communications if low cell numbers
cellchat <- filterCommunication(cellchat, min.cells = 50)

cellchat <- aggregateNet(cellchat)
cellchat@net$count
cellchat@net$weight
# Save the CellChat object
saveRDS(cellchat, file = paste0(rds@project.name, "-CC.rds"))

available_signaling <- unique(cellchat@netP$pathways)
write.csv(available_signaling, file = paste0(rds@project.name, "-CC-signals.csv"))

# Visualize interaction networks
groupSize <- as.numeric(table(cellchat@idents))
idents <- unique(cellchat@idents)

# Print network circles
pdf(paste0(rds@project.name, "-CC.pdf"), width = 8, height = 8)
par(mfrow = c(1,1), xpd=T)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 vertex.size = 4, vertex.label.cex = 1, 
                 title.name = paste(rds@project.name, "Number of Interactions"))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 vertex.size = 4, vertex.label.cex = 1,
                 title.name = paste(rds@project.name, "Interaction weights/strength"))
dev.off()

# Print separately
pdf(paste0(rds@project.name, "-CC.pdf"), width = 10, height = 6)
idents <- unique(cellchat@idents)
paths <- cellchat@netP$pathways
par(mfrow = c(2, 4), xpd=T, mar = c(1, 1, 1, 1))  # mar = c(bottom, left, top, right)
for (i in paths){
  netVisual_circle(cellchat@net$count, idents.use = c(i), vertex.weight = groupSize, 
                   vertex.size = 2, vertex.label.cex = .5, edge.width = 0.1, 
                   title.name = paste("Number of Interactions", i))
  
  netVisual_circle(cellchat@net$weight, idents.use = c(i), vertex.weight = groupSize, 
                   vertex.size = 2, vertex.label.cex = .5, weight.scale = T, 
                   title.name = paste("Interaction strength", i))
}
dev.off()

## Heatmap ####
pdf(paste0(rds@project.name, "-CCheat-",group, ".pdf"), width = 10, height = 10)
netVisual_heatmap(cellchat)
dev.off()

# Specific pathway
paths <- c("CD34", "SELE", "IL1", "PARs", "CD48", "CD200", "CEACAM", "TWEAK", "PVR",
           "ICOS", "LIGHT", "CD226", "TIGIT", "NGF", "CD137", "ICAM", "PTN", "TNF", "NPR")

pdf(paste0(cellGroup, "-CCpath-",group, ".pdf"), width = 7, height = 7)
for (path in paths){
  print(path)
  print(netVisual_heatmap(cellchat, signaling = path, color.heatmap = "Reds"))
}
dev.off()

## Circle ####
signals <- c("PECAM1", "IL16", "CLEC", "ITGB2", "TGFb")
for (signal in signals){
  png(paste0("CellChat path ", signal, ".png"))
  # Visualize signaling pathways
  netVisual_aggregate(cellchat, signaling = signal, layout = "circle")
  title(main = signal)
}
dev.off()

## Chord Diagram ####
pdf(paste0(rds@project.name, "-CCchord.pdf"), width = 30, height = 30)
# choose grid size based on number of pathways
n <- length(editPath)
nrow <- ceiling(sqrt(n))
ncol <- ceiling(n / nrow)

par(mfrow = c(nrow, ncol),
    mar = c(1, 1, 3, 1))  # small margins

for (path in editPath) {
  netVisual_aggregate(cellchat, signaling = path, layout = "chord")
  title(main = cond, cex.main = 1.2)
}
dev.off()

## Print specific signal pathways ####
signals <- c("PECAM1", "IL16", "CLEC", "ITGB2", "TGFb")
for (signal in signals){
  png(paste0("CellChat path ", signal, ".png"))
  # Visualize signaling pathways
  netVisual_aggregate(cellchat, signaling = signal, layout = "chord")
  title(main = signal)
}
dev.off()

## Bubbleplot ####
group <- "anno.geno"
idents <- levels(cellchat@idents)
paths <- cellchat@netP$pathways

pdf(paste0(rds@project.name, "-CCdot.pdf"), width = 35, height = 65)
netVisual_bubble(cellchat, #sources.use = brainNeuts, targets.use = idents, 
                 signaling = paths, 
                 title = "Brain Neutrophils Significant Pathways",
                 remove.isolate = FALSE) + coord_flip()
dev.off()

# Print by pathway
cellType <- "1-Microglia"
paths <- cellchat@netP$pathways
paths
pdf(paste0(rds@project.name, "-",group, "-CellChat-path.pdf"), width = 10, height = 6)
par(mfrow = c(2, 4), xpd=T, mar = c(1, 1, 1, 1))  # mar = c(bottom, left, top, right)
for (i in paths){
  netVisual_aggregate(cellchat, signaling = i, vertex.weight = groupSize, weight.scale = T)
}
dev.off()

## Compare Cellchat ####
c1 <- readRDS("~/nomid/Nomid-WNN-adtNorm-CC-Brain NM.rds")
c2 <- readRDS("~/nomid/Nomid-WNN-adtNorm-CC-Brain WT.rds")
c3 <- readRDS("~/nomid/Nomid-WNN-adtNorm-CC-Spleen NM.rds")
c4 <- readRDS("~/nomid/Nomid-WNN-adtNorm-CC-Spleen WT.rds")

# list(A, B) RED = UPREG in B, BLUE = DOWNREG in A
list <- list(BrainWT = c2, BrainNM= c1)
#list <- list(SpleenWT = c4, SpleenNM= c3)

cellchat <- mergeCellChat(list, add.names = names(list))

pdf("BrainNMvsWT-CC-barplot")
compareInteractions(cellchat, show.legend = F)
compareInteractions(cellchat, show.legend = F)
par(mfrow = c(1,2), xpd=TRUE)
dev.off()

pdf("BrainNMvsWT-CC.pdf", width = 10, height = 5)
netVisual_diffInteraction(cellchat, weight.scale = T, title.name = "Number of Interactions")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", title.name = "Weight of Interactions")
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
dev.off()


weight.max <- getMaxWeight(list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(list)) {
  netVisual_circle(list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(list)[i]))
}
# Monocle3 ##########
# Need to use LJI monocle3 environment for specific package downgrades
library(Seurat)
if(!require("SeuratWrappers", quietly = T)) remotes::install_github("satijalab/seurat-wrappers")
library(SeuratWrappers)
if(!require(monocle3, quietly = T)) BiocManager::install("monocle3")
library(monocle3)

# Convert Seurat object to Monocle 3 CellDataSet (CDS)
cds <- as.cell_data_set(rds)

# specfiy which reduction
reducedDims(cds)$UMAP <- Embeddings(rds, "rna.umap")

# Perform dimensionality reduction
# cds <- preprocess_cds(cds, num_dim = 50)
# cds <- reduce_dimension(cds)


# Cluster cells
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Learn trajectory graph
cds <- learn_graph(cds, use_partition = T)

# Order cells in pseudotime
# Pop-up window for root cell selection
cds <- order_cells(cds)

# Plot trajectory
pdf(paste0(rds@project.name, "CellTrajectory.pdf"))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE)
dev.off()

# DESeq ####
rds <- AddMetaData(rds, metadata = rds$patient)

pseudo <- AggregateExpression(rds, assays = "RNA", return.seurat = F, normalization.method = "none",
                              group.by = c("annotations","condition"))
counts <- pseudo$RNA
anno <- levels(unique(rds$annotations))
colData <- data.frame(sample = colnames(counts), 
                      condition = c("Act", "Act-bf", "Hlt", "Rec"),
                      annotation = anno)
rownames(colData) <- colnames(counts)

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = colData,
  design = ~ annotation + condition
)

dds <- DESeq(dds)
out <- results(dds)
tail(Cells(pseudo))
table(Idents(pseudo))
pseudo$anno.condition <- paste(pseudo$annotations, pseudo$condititon, sep = "_")

# Slingshot ####
set.seed(99)
# Slingshot ####
if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
if (!require("slingshot", quietly = TRUE)) BiocManager::install("slingshot")
library(slingshot)
if (!require("SeuratObject", quietly = TRUE)) install.packages("SeuratObject")
library(SeuratObject)
if (!require("SeuratDisk", quietly = TRUE)) remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

sce <- as.SingleCellExperiment(rds)
sce <- slingshot(sce, clusterLabels = "annotations", reducedDim = "wnn.umap")

dimred <- rds@reductions$umap@cell.embeddings
clustering <- rds$annotations
counts <- as.matrix(rds@assays$RNA@counts[rds@assays$RNA@var.features, ])
lineages <- getLineages(data = dimred, clusterLabels = clustering)

lineages