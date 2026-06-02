# CellChat ######################
# Tutorial: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#part-i-data-input-processing-and-initialization-of-cellchat-object
# Manual: https://www.rdocumentation.org/packages/CellChat/versions/1.0.0
if(!require("ComplexHeatmap", quietly = T)) BiocManager::install("ComplexHeatmap")
if (!require("CellChat", quietly = TRUE)) {
  BiocManager::install("BiocNeighbors")
  devtools::install_github("jinworks/CellChat")
}
library(CellChat)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(grid)

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
  
  # Print network circles
  pdf(paste0(sub@project.name, "-CC",cond, ".pdf"), width = 8, height = 8)
  par(mfrow = c(1,1), xpd=T)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                   vertex.size = 4, vertex.label.cex = 1, arrow.size = 0.8,
                   title.name = paste(sub@project.name, "Number of Interactions"))
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                   vertex.size = 4, vertex.label.cex = 1, arrow.size = 0.8,
                   title.name = paste(sub@project.name, "Interaction weights/strength:"))
  dev.off()
  
  # Heatmap Print by path
  paths <- cellchat@netP$pathways
  pdf(paste0(cond, "-CCheatPath.pdf"), width = 5, height = 5)
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
}

## Heatmap ####
paths <- cellchat@netP$pathways
pdf(paste0(sub@project.name, "-CCheat.pdf"), width = 5, height = 5)
netVisual_heatmap(cellchat,  title.name = paste(rds@project.name, "Number of Interactions"))
dev.off()


## Circle ####
pdf(paste0(rds@project.name, "-CCsignal.pdf"))
for (p in paths){
  netVisual_aggregate(cellchat, signaling = p, layout = "circle")
  title(paste(rds@project.name, p))
}
dev.off()

### Print by Cell Type ####
idents <- levels(cellchat@idents)
idents
pdf(paste0(rds@project.name, "-CC-cellType.pdf"), width = 8, height = 8)
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
paths <- cellchat@netP$pathways

pdf(paste0(sub@project.name, "-CCdot-",cond, ".pdf"), width = 12, height = 12)
netVisual_bubble(cellchat, #sources.use = brainNeuts, targets.use = idents, 
                 signaling = paths, title = paste(cond, "L-R Interactions"), remove.isolate = FALSE) +
  theme(axis.text.y = element_text(hjust = 0))
dev.off()


## Chord Diagram ####
### Print by path ####
pdf(paste0(rds@project.name, "-CCchord-path.pdf"), width = 16, height = 16)
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
  mtext(rds@project.name, outer = T, cex = 3)
}
dev.off()

## Compare Cellchat ####
# https://rdrr.io/github/sqjin/CellChat/f/tutorial/Comparison_analysis_of_multiple_datasets.Rmd
act <- readRDS("~/rime/CellChat/RIME-anno-CC-Act.rds")
hlt <- readRDS("~/rime/CellChat/RIME-anno-CC-Hlt.rds")
rec <- readRDS("~/rime/CellChat/RIME-anno-CC-Rec.rds")

# list(A, B) RED = UPREG in B, BLUE = DOWNREG in A
list <- list(Healthy = hlt, Active = act)

cellchat <- mergeCellChat(list, add.names = names(list))

pdf(paste0(rds@project.name, "-CC-barplot-ActvsHlt.pdf"))
compareInteractions(cellchat, show.legend = F)
compareInteractions(cellchat, show.legend = F)
par(mfrow = c(1,2), xpd=TRUE)
dev.off()

pdf(paste0(rds@project.name, "-CC-ActvsHlt.pdf"), width = 10, height = 5)
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


### Compare paths ####
actPaths <- act@netP$pathways
hltPaths <- hlt@netP$pathways
recPaths <- rec@netP$pathways

setdiff(actPaths, hltPaths)
setdiff(actPaths, hltPaths)

### ID signal groups based on functional similarity ####
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
# Requires python thru reticulate, run .rs.restartR() if new install
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

### Compare dotplot ####
netVisual_bubble(cellchat, sources.use = 4, comparison = c(1, 2), angle.x = 45)

pdf("RIME-CCdiff-ActvsHlt.pdf", width = 30, height = 10)
# Increased in set 2
netVisual_bubble(cellchat, sources.use = c(1:13), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in RIME", angle.x = 45, remove.isolate = T)
# Decreased in set 2
netVisual_bubble(cellchat, sources.use = c(1:13),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in RIME", angle.x = 45, remove.isolate = T)
dev.off()

