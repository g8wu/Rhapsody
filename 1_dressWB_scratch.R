# DRESS WB ####

## PATIENT/COND META
rds$patient <- case_when(
  rds$orig.ident == "C36-exact-poly10-51x71" ~ "DRS01",
  rds$orig.ident == "C47-1-exact-poly10-51x71" ~ "DRS01",
  rds$orig.ident == "C47-2-exact-poly10-51x71" ~ "DRS01",
  rds$orig.ident == "C85-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag03_hs" ~ "DRS02",
  rds$orig.ident == "DRS02-rec" ~ "DRS02",
  rds$orig.ident == "C86-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag05_hs" ~ "DRS03",
  rds$orig.ident == "C92-exact-poly10-51x71" ~ "DRS03",
  rds$orig.ident == "C86-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag04_hs" ~ "DRS03-BF",
  rds$orig.ident == "DRS04-act" ~ "DRS04",
  rds$orig.ident == "DRS04-rec" ~ "DRS04",
  rds$orig.ident == "C24-allfiles-exact-poly10-51x71" ~ "H01",
  rds$orig.ident == "H04" ~ "H04",
  rds$orig.ident == "C103-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag08_hs" ~ "H05",
  rds$orig.ident == "C103-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag09_hs" ~ "H06",  
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

rds$condition <- case_when(
  rds$orig.ident == "C36-exact-poly10-51x71" ~ "act",
  rds$orig.ident == "C47-1-exact-poly10-51x71" ~ "rec",
  rds$orig.ident == "C47-2-exact-poly10-51x71" ~ "rec",
  rds$orig.ident == "C85-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag03_hs" ~ "act",
  rds$orig.ident == "DRS02-rec" ~ "rec",
  rds$orig.ident == "C86-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag05_hs" ~ "act",
  rds$orig.ident == "C92-exact-poly10-51x71" ~ "rec",
  rds$orig.ident == "C86-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag04_hs" ~ "act-BF",
  rds$orig.ident == "DRS04-act" ~ "act",
  rds$orig.ident == "DRS04-rec" ~ "rec",
  rds$orig.ident == "C24-allfiles-exact-poly10-51x71" ~ "healthy",
  rds$orig.ident == "H04" ~ "healthy",
  rds$orig.ident == "C103-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag08_hs" ~ "healthy",
  rds$orig.ident == "C103-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag09_hs" ~ "healthy",  
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

table(rds$patient)
table(rds$condition)

rds$anno.cond <- paste(rds$annotations, rds$condition)
Idents(rds) <- rds$anno.cond
Idents(rds) <- factor(Idents(rds), levels = sort(levels(rds)))
rds$anno.cond <- Idents(rds)
write.csv(table(rds$patient, rds$condition),  paste0(rds@project.name, "-pat.cond.csv"))
table(rds$pat.cond)

# Mixup fix ####
Idents(rds) <- rds$patient
grep("RPS4Y", rownames(rds), value = T)
VlnPlot(rds, pt.size = 0, features = "RPS4Y1")


#  Annotate ####################
Idents(rds) <- rds$seurat_clusters
anno <- read.csv(paste0("dressWB-BASICanno-0.3.csv"), header = FALSE)
annotations <- setNames(anno[, 2], anno[, 1])
rds <- RenameIdents(rds, annotations)
rds$basicAnno <- Idents(rds)
table(Idents(rds))

# alphabetize the cell types
Idents(rds) <- factor(Idents(rds), levels = sort(levels(rds)))
rds$basicAnno <- Idents(rds)

# Barplots by Patient ####
var1 <- "annotations"
var2 <-"condition"

# patients <- levels(rds$patient)
# patients
data <- data.frame(table(rds[[var1]][,1], rds[[var2]][,1]))
colnames(data) <- c(var1, var2, "Freq")
data$Percentage <- round((data$Freq / sum(data$Freq)) * 100, 1)

# Stacked + percent
pdf(paste0(rds@project.name,"-BarPct.pdf"))
print(ggplot(data, aes(fill=condition, y=Percentage, x=annotations)) +   #EDIT
        geom_bar(position="fill", stat="identity") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Cell Composition"))
dev.off()


# Subset noUndt ####
#no small clusters 15-17
Idents(rds) <- "condition"
idents <- levels(Idents(rds))
idents
keep <- setdiff(idents, "Undetermined")
keep
rds <- subset(rds, idents = keep)
saveRDS(rds, paste0(rds@project.name,"-noUndt.rds"))

# DotPlot ####
listName <- "IL10"
group <- "anno.cond"
custom <- c(3.5, 7.5, 11.5, 15.5, 19.5, 23.5, 27.5, 31.5, 35.5)
Idents(rds) <- group
DefaultAssay(rds) <- "SCT"
genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)

# PRINT
pdf(paste(rds@project.name, listName, group, "Dot.pdf", sep = "-"), width = 11, height = 12)
for (col in names){
  print(col)
  print(DotPlot(rds, features = c("IL10"  ,    "IL10RA" ,   "IL10RB",    "IL10RB-DT"), cols = "RdYlBu", col.min = 0, dot.scale = 5) +
          #coord_flip() + 
          # ggtitle(col) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) +
          # geom_hline(yintercept = seq(n-1.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black", linetype = "dashed") +
          geom_hline(yintercept = custom, color = "black")
        # geom_hline(yintercept = custom, color = "black")
  )
}
dev.off()

# CellChat ####
library(CellChat)
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
library(patchwork)
setwd("~/dress/spatial/CellChat/500_cell_cutoff/")
healthy <- readRDS("~/dress/spatial/CellChat/500_cell_cutoff/Dress_healthy_BeC-CC-Healthy.rds")
dress <- readRDS("~/dress/spatial/CellChat/500_cell_cutoff/Dress_healthy_BeC-CC-DRESS.rds")

# list(A, B) RED = UPREG in B, BLUE = DOWNREG in A
list <- list(Healthy = healthy, DRESS= dress)
list <- list(Healthy = healthy, DRESS_active= act)
list <- list(Healthy = healthy, DRESS_recovered= rec)
list <- list(DRESS_active = act, DRESS_BF= BF)

cellchat <- mergeCellChat(list, add.names = names(list))

compareInteractions(cellchat, show.legend = F)
compareInteractions(cellchat, show.legend = F)
par(mfrow = c(1,2), xpd=TRUE)

pdf("Dress-CC-diff-ActvBF.pdf", width = 10, height = 5)
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


# Abseq specific ####
monoKeep <- c("CD14", "CD16-FCGR3A", "CD11b", "CD33", "CD64-FCGR1A", "HLA-DR-CD74",
              "CD32-FGCR2A", "CD11c", "CD44", "CD62L", "CD86",  "CD184-CXCR4")
rds <- readRDS("~/dress/WB/subclusters_adtNorm/Monocyte/Monocyte-RNAADT.rds")
adt.matrix <- rds[["ADT"]]$counts
new.adt.matrix <-  adt.matrix[monoKeep,]
rds[["ADT"]] <- CreateAssayObject(counts= new.adt.matrix)

# Set adtNorm positive peak to CD33
set.seed(99)
# Make new folder for each subcluster
print(cluster)
dir.create(paste0(wkdir, "/", cluster))
setwd(paste0(wkdir,"/", cluster))

## ADTnorm #####
print("____________________________________________ ADTnorm start")
cell_x_adt <- t(rds@assays$ADT@counts)
cell_x_feature <- data.frame(rds@active.ident)

# Make sure there is a 'sample' column set to an ident
colnames(cell_x_feature) <- 'sample'
head(cell_x_feature)

cell_x_adt_norm = ADTnorm(
  cell_x_adt = cell_x_adt, 
  cell_x_feature = cell_x_feature, 
  save_outpath = getwd(), 
  study_name = rds@project.name,
  marker_to_process = rownames(rds@assays$ADT),
  bimodal_marker = NULL,             # default NULL: try different settings to find bimodal peaks
  # trimodal_marker = c("CD45-F11-Ptprc"),       # CD4 and CD45RA tend to have 3 peaks
  # setting the CD3 uni-peak of buus_2021_T study to positive peak if only one peak is detected for CD3 marker
  # positive_peak = list(ADT = "CD3", sample = "buus_2021_T"), 
  positive_peak = list(ADT = "CD3"), 
  brewer_palettes = "Dark2",
  save_fig = TRUE,
  target_landmark_location = "fixed",
  shoulder_valley = T,               # Look for "shoulder" as pos peak (technical variation -> no clear separation b/w neg/pos)
  #multi_sample_per_batch = T,        # Omit aligning the one pos peak
  #customize_landmark = T             # Manual adjustment UI 
)

# Put ADTNorm matrix back into rds
rds@assays$ADT@scale.data <- t(cell_x_adt_norm)
print(paste(Sys.time(), "ADTnorm done"))

## PCA ####
rds <- RunPCA(rds, assay = "SCT", reduction.name = "pca_rna")
VariableFeatures(rds@assays$ADT) <- rownames(rds@assays$ADT)
rds <- RunPCA(rds, assay = "ADT", reduction.name = "pca_adt", npcs = adtPCs)

## Harmonize RNA ####
# Batch effect correct (BeC)
rds <- RunHarmony(rds, group.by.vars = "orig.ident", 
                  reduction = "pca_rna", 
                  reduction.save = "harmony_rna")
### UMAP RNA ####
# UMAP Without Harmony
rds <- RunUMAP(rds, reduction = "pca_rna", dims = 1:rnaPCs, reduction.name = "umap")
# UMAP With Harmony
rds <- RunUMAP(rds, reduction = "harmony_rna", dims = 1:rnaPCs, reduction.name = "harmony_umap")

# Print
pdf(paste0(rds@project.name,"-PCA-Harmony.pdf"))
Idents(rds) <- "orig.ident"
print(DimPlot(rds, reduction = "umap"))
print(DimPlot(rds, reduction = "harmony_umap"))
print(DimPlot(rds, reduction = "pca_rna"))
print(DimPlot(rds, reduction = "pca_adt"))
dev.off()

## Multi-res Clust ####
### RNA ONLY ####
print("_____________________ Starting RNA ONLY Clustering")
set.seed(99)
rds <- FindNeighbors(rds, dims = 1:rnaPCs, reduction = "harmony_rna")

# Multi-res Clustering
pdf(paste0(rds@project.name, "-UMAP-RNA.pdf"))
for(res in resRange) {
  print(res)
  rds <- FindClusters(rds, algorithm = alg, resolution = res)
  print(DimPlot(rds, reduction = 'harmony_umap', label = TRUE, raster = F) +
          labs(title = paste0("RNA UMAP ", algKey[alg],": ", res)))
}
dev.off()

print(paste0(Sys.time(), " -> RNA ONLY clustering done!"))

### RNA + ADT ####
print("_____________________ Starting RNA & ADT Clustering")
set.seed(99)
rds <- FindMultiModalNeighbors(rds, reduction.list = list("harmony_rna", "pca_adt"),
                               dims.list = list(1:rnaPCs, 1:adtPCs))
rds <- RunUMAP(rds, nn.name= "weighted.nn", reduction.name = "wnn.umap")
# Multi-res Clustering
pdf(paste0(rds@project.name, "-UMAP-RNAADT.pdf"))
for(res in resRange) {
  print(res)
  rds <- FindClusters(rds, algorithm = alg, resolution = res, graph.name = "wsnn")
  print(DimPlot(rds, reduction= "wnn.umap", label = T, raster = F) + 
          labs(title = paste0("ADT&RNA UMAP ", algKey[alg],": ", res)))
}
dev.off()

## SAVE!! ####
saveRDS(rds, file= paste0(rds@project.name, "-RNAADT.rds"))

print(paste0(Sys.time(), " -> RNA & ADT clustering done!"))

## QC FeatPlot ####
features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "ADT_total")
pdf(paste0(rds@project.name, "-QCUMAP.pdf"))
print(FeaturePlot(rds, reduction = "harmony_umap", features = features, ncol = 2))
print(FeaturePlot(rds, reduction = "wnn.umap", features = features, ncol = 2))
dev.off()

## Umap orig.ident ####
pdf(paste(rds@project.name, "UMAP-batchEffectCheck.pdf", sep = "-"))
DimPlot(rds, reduction = "harmony_umap", group.by = "orig.ident")
DimPlot(rds, reduction = "harmony_umap", group.by = "pat.cond")
DimPlot(rds, reduction = "wnn.umap", group.by = "orig.ident")
DimPlot(rds, reduction = "wnn.umap", group.by = "pat.cond")
dev.off()

## Abseq Featplot #### 
abseq <- rownames(rds@assays$ADT)

# PRINT
pdf(paste(rds@project.name, "AbsFeats-RNA.pdf", sep = "-"), width = length(abseq) / 5 * 4.8, height = length(abseq) / 5 * 5)
rds@misc$umap <- "harmony_umap"
DefaultAssay(rds) <-"ADT"
plots <- FeaturePlot(rds, reduction = rds@misc$umap, features = abseq, ncol = 5) & 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank())
DefaultAssay(rds) <-"SCT"
print(plots+ plot_annotation(title ="Abseq", theme = theme(plot.title = element_text(size =50))))
dev.off()
pdf(paste(rds@project.name, "AbsFeats-WNN.pdf", sep = "-"), width = length(abseq) / 5 * 4.8, height = length(abseq) / 5 * 5)
rds@misc$umap <- "wnn.umap"
DefaultAssay(rds) <-"ADT"
plots <- FeaturePlot(rds, reduction = rds@misc$umap, features = abseq, ncol = 5) & 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank())
DefaultAssay(rds) <-"SCT"
print(plots+ plot_annotation(title ="Abseq", theme = theme(plot.title = element_text(size =50))))
dev.off()

# Subcluster pt2 ####
rds@misc$umap <- "wnn.umap"
rds$seurat_clusters <- rds$wsnn_res.0.05
rds@project.name <-"Monocyte-testRmvADT-WNN"

# rds@misc$umap <- "harmony_umap"
# rds@project.name <-"Neutrophil-RNA"
# rds$seurat_clusters <- rds$SCT_snn_res.0.1

## Subset umap ####
group <- "seurat_clusters"
split.by <- "condition"

pdf(paste(rds@project.name, "UMAPsubset.pdf", sep = "-"))
DimPlot(rds, reduction = rds@misc$umap, raster = F, group.by = group,split.by = split.by, ncol=2) + 
  NoLegend()
dev.off()



## VlnPlot ####
group <- "seurat_clusters"
Idents(rds) <- group

pdf(paste0(rds@project.name, "-", group, "-QCVln", ".pdf"), height = 15, width = 10)
plots <- lapply(features[1:3], function(i) {
  VlnPlot(rds, features = i, pt.size = 0)
})
# Switch to ADT to print ADT_total
adtPlot <- VlnPlot(rds, features = features[4], pt.size = 0)
plots <- c(plots, list(adtPlot))
print(wrap_plots(plots, ncol = 2) & NoLegend() & 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)))
dev.off()

## Abseq Ridge ####
Idents(rds) <- rds$seurat_clusters
DefaultAssay(rds) <- "ADT"
abseq <- rownames(rds@assays$ADT)
plots <- RidgePlot(rds, features = abseq, stack = T) + NoLegend()
# PRINT
pdf(paste0(rds@project.name, "-AbseqRidge-", group, ".pdf"))
print(plots + labs(title = paste("Abseq Ridge Plot:", group)) + 
        theme(axis.text.y = element_text(hjust = 0))
)
dev.off()
DefaultAssay(rds) <- "SCT"

## Abseq Feat ####
rds$clust.cond <- paste(rds$seurat_clusters, rds$condition)

clustSize <- table(Idents(rds))
keep  <- names(clustSize[clustSize >= 20])
sub <-subset(rds, idents = keep)
group <- "clust.cond"
abseq <- rownames(rds[["ADT"]])
DefaultAssay(sub) <- "ADT"
pdf(paste(rds@project.name, group, "AbseqDot.pdf", sep = "-"),width = 10, height = 12)
DotPlot(sub, features = abseq, col.min = 0, cols = "RdYlBu", group.by = group) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(hjust = 0)) + 
  ggtitle("AbSeq") +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5)
dev.off()
DefaultAssay(rds) <- "SCT"

## Abseq Dot ####
group <-"seurat_clusters"
Idents(rds) <- group
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
pdf(paste(rds@project.name, group, "AbseqDot.pdf", sep = "-"),width = 10, height = 6)
DotPlot(rds, features = abseq, col.min = 0, cols = "RdYlBu", group.by = group) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(hjust = 0)) + 
  ggtitle("AbSeq") +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5)
dev.off()
DefaultAssay(rds) <- "SCT"

## Barplots ####
table <- data.frame(table(rds$seurat_clusters, rds$condition))
table$Percentage <- round((table$Freq / sum(table$Freq)) * 100, 1)

pdf(paste0(rds@project.name, "-Bar.pdf"), width = 5, height = 5)
ggplot(table, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_text(aes(label = Freq), # Show count % Pct
  #           vjust = -0.5, size =1) +
  labs(title = "Cell Type Distribution", x = "Cluster", y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
dev.off()

var1 <- "seurat_clusters"
var2 <-"condition"
data <- data.frame(table(rds[[var1]][,1], rds[[var2]][,1]))
colnames(data) <- c(var1, var2, "Freq")
data$Percentage <- round((data$Freq / sum(data$Freq)) * 100, 1)

# Write to csv file
write.csv(table(rds[[var1]][,1], rds[[var2]][,1], rds$patient), paste0(rds@project.name,"-table.csv"))

# Stacked + percent
pdf(paste0(rds@project.name,"-BarPct.pdf"))
print(ggplot(data, aes(fill=condition, y=Percentage, x=seurat_clusters)) +   #EDIT
        geom_bar(position="fill", stat="identity") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Cell Composition"))
dev.off()

## DEG ####
### Subcluster ####
group <- "seurat_clusters"
Idents(rds)<-group
idents <- levels(Idents(rds))
idents
degs <- FindAllMarkers(rds, only.pos = F, logfc.threshold = 0.5, assay = "SCT", recorrect_umi = F)
dim(degs)

# Create Excel workbook
wb <- createWorkbook()
for (i in idents){
  print(i)
  subset <- degs[degs$cluster == i,]
  # Sort by decreasing avg LogFC
  subset <- subset[order(subset$avg_log2FC, decreasing = T),]
  addWorksheet(wb, i)
  writeData(wb, sheet = i , subset)
}
saveWorkbook(wb, file=paste0(rds@project.name, "-", group, "-DEGs.xlsx"), overwrite = TRUE)

### Conditions ####
degs <- NULL
wb <-NULL
group <- "condition"
Idents(rds) <- group
comparisons <- list(list("act","rec"), list("act", "healthy"), list("rec", "healthy"),
                    list("act-BF", "act"), list("act-BF", "healthy"), list("act-BF", "rec"))
#rds <- PrepSCTFindMarkers(rds)
wb <- createWorkbook()
for (pair in comparisons){
  A = pair[[1]]
  B = pair[[2]]
  degs <- FindMarkers(rds, only.pos = F, ident.1 = A, ident.2 = B,
                      min.pct = 0.1, logfc.threshold = 0.15, assay = "SCT")
  degs$genes <- rownames(degs)
  degs <- degs[order(degs$avg_log2FC, decreasing = T),]
  addWorksheet(wb, sheet = paste0(A,".", B))
  writeData(wb, sheet = paste0(A,".", B), degs)
}
saveWorkbook(wb,file=paste0(rds@project.name, "-", group, "-DEGs.xlsx"), overwrite = TRUE)

#### Filter DEGs ####
# by LogFC, gate out upregs if pct1 < 0.5
## IF BLANK SHEETS MANUALLY DELETE!!!!!! ####
input_file_list <- c(paste0(rds@project.name, "-seurat_clusters-DEGs"), paste0(rds@project.name, "-condition-DEGs"))

for (input_file in input_file_list){
  print(input_file)
  sheet_names <- excel_sheets(paste0(input_file, ".xlsx"))
  
  wb <- createWorkbook()
  for (sheet in sheet_names) {
    # Read the sheet
    data <- read_excel(paste0(input_file, ".xlsx"), sheet = sheet)
    
    # Apply filtering logic
    top <- data %>%
      filter(p_val_adj < 0.05 & (avg_log2FC > 0 & pct.1 >= 0.5)) %>% 
      slice_head(n=100)
    bot <- data %>%
      filter(p_val_adj < 0.05 & (avg_log2FC < 0 & pct.2 >= 0.5)) %>% 
      slice_tail(n=100)
    
    all <- bind_rows(top, bot)
    # Add to workbook
    addWorksheet(wb, sheetName = sheet)
    writeData(wb, sheet = sheet, all)
  }
  saveWorkbook(wb, file=paste0(input_file, "-filtered.xlsx"), overwrite = TRUE)
  
  #### Dotplot filtered DEGs ####
  print("Dotplot")
  colors <- colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))
  pdf(paste(input_file, "DotTop100.pdf", sep = "-"), width = 20, height = 6)
  for (sheet in sheet_names) {
    # Read the sheet
    degs <- read.xlsx(paste0(input_file, "-filtered.xlsx"), sheet = sheet)
    # Print cluster name
    print(sheet)
    plot.new()
    # x and y coordinates are from 0 (left/bottom) to 1 (right/top)
    text(x = 0.5, y = 0.5, paste(sheet, "100 Upreg"), cex = 5, font = 2) # cex changes text size, font changes style
    
    Clustered_DotPlot(rds, features = degs$gene, flip = T, 
                      colors_use_exp = rev(colors(20)), x_lab_rotate = T, 
                      group.by = "seurat_clusters", exp_color_min = 0, 
                      plot_km_elbow = F)
  }
  dev.off()
  
  #### DEG.xlsx -> Volcano ####
  print("Volcano")
  pdf(paste0(input_file, "-Volcano.pdf"))
  for (sheet in sheet_names){
    print(sheet)
    degs <- read.xlsx(paste0(input_file, ".xlsx"), sheet = sheet)
    print(EnhancedVolcano(degs, lab = degs$gene, x = "avg_log2FC", y = "p_val_adj",
                          pCutoff = 0.05, FCcutoff = 2.0,
                          title = sheet, subtitle = ""))
  }
  dev.off()
}


## Dotplot top 10 DEGs ####
# FILTER OUT SMALL CLUSTERS
n <- 50
rds$clust.cond <- paste(rds$seurat_clusters, rds$condition)
Idents(rds) <- rds$clust.cond
clustSize <- table(rds$clust.cond)
keep  <- names(clustSize[clustSize >= n])
sub <-subset(rds, idents = keep)
#### by cluster ####
input_file <- paste0(paste0(rds@project.name, "-seurat_clusters-DEGs"))
sheet_names <- excel_sheets(paste0(input_file, "-filtered.xlsx"))
genes <- c()
for (sheet in sheet_names) {
  print(sheet)
  # Read the sheet
  degs <- read.xlsx(paste0(input_file, "-filtered.xlsx"), sheet = sheet)
  genes <- append(genes, degs$gene[1:10])
}
genes <- unique(genes[!is.na(genes)])
genes <- genes[!grepl("^ENS", genes)]
genes <- genes[!grepl("^LINC", genes)]
pdf(paste0(rds@project.name, "-seurat_clusters-DotTop10.pdf"))
print(DotPlot(sub, features = unique(genes), cols = "RdYlBu", col.min = 0, dot.scale = 5, group.by = "clust.cond") +
        ggtitle(paste(rds@project.name, "Top 10 per cluster: clusters smaller than",n, "cells filtered out")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) +
        geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5))
dev.off()

#### by condition ####
# FILTER OUT SMALL CLUSTERS
n <- 50
rds$cond.clust <- paste(rds$condition, rds$seurat_clusters)
Idents(rds) <- "cond.clust"
clustSize <- table(rds$cond.clust)
keep  <- names(clustSize[clustSize >= n])
sub <-subset(rds, idents = keep)

input_file <-paste0(rds@project.name, "-condition-DEGs")
sheet_names <- excel_sheets(paste0(input_file, ".xlsx"))
genes <- c()
for (sheet in sheet_names) {
  print(sheet)
  # Read the sheet
  degs <- read.xlsx(paste0(input_file, "-filtered.xlsx"), sheet = sheet)
  genes <- append(genes, degs$gene[1:10])
}
genes <- unique(genes[!is.na(genes)])
genes <- genes[!grepl("^ENS", genes)]
genes <- genes[!grepl("^LINC", genes)]

pdf(paste(input_file, "DotTop10.pdf", sep = "-"), width = 10, height = 7)
print(DotPlot(sub, features = unique(genes), cols = "RdYlBu", col.min = 0, dot.scale = 5, group.by = "cond.clust") +
        ggtitle(paste(rds@project.name, "Top 10 per genotype: clusters smaller than",n, "cells filtered out")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) +
        geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5))
dev.off()

