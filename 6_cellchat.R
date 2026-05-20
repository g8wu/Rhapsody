# CellChat ####
if (!require("CellChat", quietly = TRUE)) devtools::install_github("jinworks/CellChat")
library(CellChat)
if(!require("ComplexHeatmap", quietly = T)) BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(grid)
library(circlize)

Idents(rds) <- "condition"
idents <- levels(Idents(rds))
idents

for (cond in idents){
  print(cond)
  sub <- subset(rds, idents = cond)
  Idents(sub) <- "annotations"
  DefaultAssay(sub) <- "SCT"
  data.input <- sub[["SCT"]]$data
  labels <- Idents(sub)
  meta <- data.frame(group = labels, row.names = names(labels))
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  
  # Set CellChat database
  #CellChatDB <- CellChatDB.mouse  # mouse data
  CellChatDB <- CellChatDB.human  # human data
  showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChatDB
  
  # Subset data
  cellchat <- subsetData(cellchat)
  
  # Identify overexpressed genes and interactions
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probability
  cellchat <- computeCommunProb(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Filter out communications if low cell numbers
  cellchat <- filterCommunication(cellchat, min.cells = 50)
  
  cellchat <- aggregateNet(cellchat)
  cellchat@net$count
  cellchat@net$weight
  paths <- unique(cellchat@netP$pathways)
  write.csv(paths, file = paste0(sub@project.name, "-",cond, "-CCsignals.csv"))
  
  # Save the CellChat object
  saveRDS(cellchat, file = paste0(sub@project.name, "-CC-",cond, ".rds"))
  
  # Visualize interaction networks
  groupSize <- as.numeric(table(cellchat@idents))
  idents <- unique(cellchat@idents)
  
  ## CirclePlots ####
  pdf(paste0(sub@project.name, "-CC-",cond, ".pdf"), width = 8, height = 8)
  par(mfrow = c(1,1), xpd=T)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                   vertex.size = 4, vertex.label.cex = 1, arrow.size = 0.8,
                   title.name = paste(cond, "Number of Interactions"))
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                   vertex.size = 4, vertex.label.cex = 1, arrow.size = 0.8,
                   title.name = paste(cond, "Interaction weights/strength:"))
  dev.off()
  
  ## Heatmap ####
  paths <- cellchat@netP$pathways
  pdf(paste0(sub@project.name, "-CCheat-",cond, ".pdf"), width = 5, height = 5)
  netVisual_heatmap(cellchat,  title.name = paste(cond, "Number of Interactions"))
  dev.off()
  
  # Print by path
  pdf(paste0(sub@project.name, "-CCheatPath-",cond, ".pdf"), width = 5, height = 5)
  for(path in paths) {
    print(path)
    tryCatch({
      heats <- netVisual_heatmap(cellchat, signaling = path, title.name = paste(cond, path, ": Number of Interactions"))
      draw(heats)  
    }, error = function(x) {
      message("Skipping: ", conditionMessage(x))
    })
  }
  dev.off()
  
  ## Circle by Signal ####
  pdf(paste0(sub@project.name, "-CCsignal-", cond,".pdf"))
  for (p in editPath){
    netVisual_aggregate(cellchat, signaling = p, layout = "circle")
    title(paste(cond, p))
  }
  dev.off()
}


### Print by Cell Type ####
idents <- levels(cellchat@idents)
idents
pdf(paste0(rds@project.name, "-CC-cellType-", cond, ".pdf"), width = 8, height = 8)
# choose grid size based on number of pathways
n <- length(idents)
nrow <- ceiling(sqrt(n))
ncol <- ceiling(n / nrow)

par(mfrow = c(nrow, ncol),
    mar = c(0, 0, 0, 0),
    oma = c(0, 0, 3, 0))  # small margins
for(ident in idents) {
  print(ident)
  tryCatch({
    netVisual_aggregate(cellchat, layout = "circle", signaling = paths, sources.use = ident)
  }, error = function(x) {
    message("Skipping: ", conditionMessage(x))
  })
  mtext(cond, outer = T, cex = 3)
}
dev.off()

## CC Dot ####
pdf(paste0(sub@project.name, "-CCdot-",cond, ".pdf"), width = 10, height = 12)
# choose grid size based on number of pathways
n <- length(editPath)
nrow <- ceiling(sqrt(n))
ncol <- ceiling(n / nrow)
par(mfrow = c(nrow, ncol), mar = c(0, 0, 0, 0), oma = c(0, 0, 3, 0))
for (path in editPath) {
  netVisual_bubble(cellchat, #sources.use = brainNeuts, targets.use = idents, 
                   signaling = paths, title = paste(cond, path, "L-R Interactions"), remove.isolate = FALSE) +
    theme(axis.text.y = element_text(hjust = 0))
}
mtext(cond, outer = T, cex = 3)
dev.off()


## Chord Diagram ####
### Print by path ####
pdf(paste0(rds@project.name, "-CCchord-path-", cond, ".pdf"), width = 16, height = 16)
# choose grid size based on number of pathways
n <- length(paths)
nrow <- ceiling(sqrt(n))
ncol <- ceiling(n / nrow)

par(mfrow = c(nrow, ncol),
    mar = c(0, 0, 0, 0),
    oma = c(0, 0, 3, 0))  # small margins
for(path in paths) {
  print(path)
  tryCatch({
    netVisual_aggregate(cellchat, signaling = path, layout = "chord" )  
  }, error = function(x) {
    message("Skipping: ", conditionMessage(x))
  })
  mtext(cond, outer = T, cex = 3)
}
dev.off()

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
