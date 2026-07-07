# RIME SCRATCH ####
set.seed(99)
library(Seurat)
library(dplyr)
n <- 50     # Let variables reach up to n GB
options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB
# Metadata ####
C39 <- readRDS("~/Projects/Rime/samples/C39-allfiles-basic-poly10-51x71_Seurat.rds")
C98 <- readRDS("~/Projects/Rime/samples/C98-exact-poly10-51x71_Seurat.rds")
C114 <- readRDS("~/Projects/Rime/samples/C114-expected-poly10-51x71_Seurat.rds")
C119 <- readRDS("~/Projects/Rime/samples/C119-exact-poly10-51x71_Seurat.rds")
C120 <- readRDS("~/Projects/Rime/samples/C120-exact-poly10-51x71_Seurat.rds")
## Patient ####
C39$patient <- case_when(
  C39$orig.ident == "C39-allfiles-basic-poly10-51x71_Seurat.rds" ~ "RIM01"
)

C98$patient <- case_when(
  C98$Sample_Tag == "SampleTag12_hs" ~ "RIM02",
  C98$Sample_Tag == "SampleTag01_hs" ~ "RIM02",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

C119$patient <- case_when(
  C119$Sample_Tag == "SampleTag06_hs" ~ "RIM04",
  C119$Sample_Tag == "SampleTag07_hs" ~ "RIM04",
  C119$Sample_Tag == "SampleTag08_hs" ~ "RIM05",
  C119$Sample_Tag == "SampleTag09_hs" ~ "RIM05",
  C119$Sample_Tag == "SampleTag10_hs" ~ "RIM05",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

C120$patient <- case_when(
  C120$Sample_Tag == "SampleTag06_hs" ~ "RIM03",
  C120$Sample_Tag == "SampleTag07_hs" ~ "RIM03",
  C120$Sample_Tag == "SampleTag08_hs" ~ "RIM06",
  C120$Sample_Tag == "SampleTag09_hs" ~ "RIM06",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

C114$patient <- case_when(
  C114$Sample_Tag == "SampleTag01_hs" ~ "H01",
  C114$Sample_Tag == "SampleTag02_hs" ~ "H04",
  C114$Sample_Tag == "SampleTag03_hs" ~ "H05",
  C114$Sample_Tag == "SampleTag04_hs" ~ "H06",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)
## Condition ####
C39$condition <- case_when(
  C39$Sample_Tag == "SampleTag01_hs" ~ "Act",
  C39$Sample_Tag == "SampleTag03_hs" ~ "Rec",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

C98$condition <- case_when(
  C98$Sample_Tag == "SampleTag12_hs" ~ "Act",
  C98$Sample_Tag == "SampleTag01_hs" ~ "Rec",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

C119$condition <- case_when(
  C119$Sample_Tag == "SampleTag06_hs" ~ "Act",
  C119$Sample_Tag == "SampleTag07_hs" ~ "Rec",
  C119$Sample_Tag == "SampleTag08_hs" ~ "Act",
  C119$Sample_Tag == "SampleTag09_hs" ~ "Rec",
  C119$Sample_Tag == "SampleTag10_hs" ~ "Act",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

C120$condition <- case_when(
  C120$Sample_Tag == "SampleTag06_hs" ~ "Act",
  C120$Sample_Tag == "SampleTag07_hs" ~ "Rec",
  C120$Sample_Tag == "SampleTag08_hs" ~ "Act",
  C120$Sample_Tag == "SampleTag09_hs" ~ "Rec",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

C114$condition <- case_when(
  C114$Sample_Tag == "SampleTag01_hs" ~ "Hlt",
  C114$Sample_Tag == "SampleTag02_hs" ~ "Hlt",
  C114$Sample_Tag == "SampleTag03_hs" ~ "Hlt",
  C114$Sample_Tag == "SampleTag04_hs" ~ "Hlt",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)


# Merge and subset ####
list <- list(C39, C98, C114, C119, C120)
rds <- merge(x = list[[1]], y = list[-1], project = "RIME")
Idents(rds) <- rds$condition
table(Idents(rds))
table(rds$condition)
table(rds$patient)

# Remove Undetermined
rds <- subset(rds, idents = c("Rec", "Act", "Hlt"))

# Sample Type ####
rds$sample_type <- case_when(
  rds$orig.ident == "C119-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_hs" ~ "BF",
  TRUE ~ "PBMC"
)
table(rds$sample_type)

# SAVE!!
saveRDS(rds, "Rime.rds")
# Select res ####
rds@misc$umap <- "harmony_umap"
rds$seurat_clusters <- rds$SCT_snn_res.0.5
Idents(rds) <- rds$seurat_clusters
# order by numeric value
Idents(rds) <- factor(Idents(object = rds), levels = sort(as.numeric(levels(rds))))
rds$seurat_clusters <- Idents(rds)

# UMAPs post annotation ####
pdf(paste0(rds@project.name,"-UMAP.pdf"))
DimPlot(rds, reduction = rds@misc$umap, label = T,raster = F, group.by = "annotations", repel = 2) + NoLegend()
DimPlot(rds, reduction = rds@misc$umap,raster = F, group.by = "annotations", split.by = "condition", ncol = 2) + NoLegend()
dev.off()

# CellChat ####
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
  cellchat@idents <- factor(cellchat@idents, levels = all.types)
  
  # Set CellChat database
  #CellChatDB <- CellChatDB.mouse  # mouse data
  CellChatDB <- CellChatDB.human  # human data
  showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChatDB
  
  # Keep all identities
  all.types <- levels(cellchat@idents)
  cellchat@idents <- factor(cellchat@idents, levels = all.types)
  cellchat@idents <- droplevels(cellchat@idents)
  
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
  pdf(paste0(sub@project.name, "-CC-",cond, ".pdf"), width = 8, height = 8)
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
# Check sources index
levels(cellchat@ident)
pdf(paste0(rds@project.name, "-CCsignal.pdf"))
for (p in paths){
  netVisual_aggregate(cellchat, signaling = p, layout = "circle")
  title(paste(rds@project.name, p))
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



# UMAPS ####
pdf(paste0("RIME-RNAonly-UMAPs.pdf"))
DimPlot(clr, reduction = "harmony_umap",label= T) + ggtitle("RNA only: res 0.2, Leiden")
DimPlot(clr, reduction = "harmony_umap", group.by = "orig.ident")
DimPlot(clr, reduction = "harmony_umap", split.by = "orig.ident", ncol = 2) +NoLegend()
DimPlot(clr, reduction = "harmony_umap", group.by = "condition")
DimPlot(clr, reduction = "harmony_umap", split.by = "condition", ncol = 2) +NoLegend()
DimPlot(clr, reduction = "harmony_umap", group.by = "patient")
DimPlot(clr, reduction = "harmony_umap", split.by = "patient", ncol = 3) +NoLegend()
dev.off()
