# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
# BiocManager::install('multtest')
set.seed(99)
library(dplyr)
if(!require("scCustomize", quietly = T)) install.packages("scCustomize")
if(!require("openxlsx", quietly = T)) install.packages("openxlsx")
if(!require("readxl", quietly = T)) install.packages("readxl")
if(!require("EnhancedVolcano", quietly = T)) BiocManager::install("EnhancedVolcano")
n
if(!require("scCustomize", quietly = T)) install.packages("scCustomize")
if(!require("ComplexHeatmap", quietly = T)) BiocManager::install("ComplexHeatmap")
n
if (!require("harmony", quietly = TRUE)) BiocManager::install("harmony")
n
BiocManager::install("glmGamPoi")
n
if (!require("ADTnorm", quietly = TRUE)) { BiocManager::install(c("cytolib","flowCore", "ncdfFlow", "flowViz", "flowWorkspace","flowStats"))
  remotes::install_github("yezhengSTAT/ADTnorm", build_vignettes = FALSE)
}
n


library(scDblFinder)
library(ggplot2)
library(sctransform)
library(purrr)
library(harmony)
library(Seurat)
library(ADTnorm)
library(glmGamPoi)
library(RColorBrewer)
#library(edgeR)
library(ggrepel)
library(gridExtra)
library(scCustomize)
library(patchwork)
library(tidyverse)
library(dplyr)
library(SingleCellExperiment)
library(readxl)
library(openxlsx)
library(EnhancedVolcano)
# if(!require("metap", quietly = T)) install.packages("metap")
# library(metap)
display.brewer.all(colorblindFriendly = TRUE)
setwd("/mnt/bioadhoc/Groups/Collaborators/ben.croker/nomid")

#  Annotate ####################
Idents(rds) <- rds$seurat_clusters
anno <- read.csv(paste0("nomid-anno.csv"), header = FALSE)
annotations <- setNames(anno[, 2], anno[, 1])
rds <- RenameIdents(rds, annotations)
rds$annotations <- Idents(rds)
table(Idents(rds))

# alphabetize the cell types
Idents(rds) <- factor(Idents(rds), levels = sort(levels(rds)))
rds$basicAnno <- Idents(rds)

# order by numeric value
Idents(rds) <- factor(Idents(object = rds), levels = sort(as.numeric(levels(rds))))

rds@project.name <- paste0(rds@project.name, "-Anno")
saveRDS(rds, paste0(rds@project.name, ".RDS"))

## Barplot #########
library(gridExtra)
# PRINT
pdf(paste0(rds@project.name, "-Bar.pdf"), width = 5, height = 5)
table <- data.frame(table(rds$annotations, rds$tiss.geno.formal))
table$Percentage <- round((table$Freq / sum(table$Freq)) * 100, 1)
ggplot(table, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_text(aes(label = Freq), # Show count % Pct
  #           vjust = -0.5, size =1) +
  labs(title = "Cell Type Distribution", x = "Cluster", y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

# Barplot Percentage####
var1 <- "annotations"
var2 <-"tiss.geno.formal"
data <- data.frame(table(rds[[var1]][,1], rds[[var2]][,1]))
colnames(data) <- c(var1, var2, "Freq")
data$Percentage <- round((data$Freq / sum(data$Freq)) * 100, 1)

# Stacked + percent
print(ggplot(data, aes(fill=tiss.geno.formal, y=Percentage, x=annotations)) +   #EDIT
        geom_bar(position="fill", stat="identity") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle("Cell Composition"))
dev.off()

### Stacked + Pct ####
ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity")
var1 <- "annotations"
var2 <-"geno"
data <- data.frame(table(rds[[var1]][,1], rds[[var2]][,1]))

ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity")

plot <- ggplot(data) + aes(fill=annotations, y=Percentage, x=geno) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~ Condition) +NoLegend()

pdf(paste0(rds@project.name,"-BarStackPct.pdf"))
print(plot)
dev.off()

### Unique UMI ####
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
# PRINT
pdf(paste(rds@project.name, "UMIBarplot.pdf", sep = "-"))
print(umiBar)
dev.off()

# Cell Select Box ####
focus <- c("Microglia")

# Get target cluster
sub <- subset(rds, idents = focus)

# WhichCells subset check
sample <- WhichCells(rds, idents = focus)
DimPlot(rds, reduction = sub@misc$umap, cells.highlight = sample, label = T, raster=F) + 
  NoLegend()

# Individual cells
# HoverLocator(DimPlot(sub, reduction = sub@misc$umap,raster = F),
#              information = FetchData(sub, vars = c("ident", "PC_1", "nFeature_RNA")))
# select.cells <- c("3112740_5")

# Box selection
sub <-rds
select.cells <- CellSelector(DimPlot(sub, reduction = rds@misc$umap,raster = F))
# write.table(select.cells, file = "Eos_Log.csv", append = TRUE, 
#             sep = ",", col.names = FALSE, row.names = FALSE)
# Assign reselection to subset 
# Idents(sub, cells = select.cells) <- "1-Microglia"

# Check reselection on subset
# DimPlot(sub, reduction = rds@misc$umap, label = T)

# Set reselection to original object
Idents(rds, cells = select.cells) <- "Endothelial"
Idents(rds) <- factor(Idents(rds), levels = sort(levels(rds)))
#Idents(rds) <- factor(Idents(rds), levels = sort((levels(rds))))
rds$basicAnno <- Idents(rds)
DimPlot(rds, reduction = rds@misc$umap, label = T)
gc()

# DEG: Subclust A vs B #######
A = "B-like Neutrophil"
B = "Pre/Immature Neutrophil"
group = "seurat_clusters"
Idents(rds) <- group

rds <- PrepSCTFindMarkers(rds)
degs <- FindMarkers(rds, only.pos = F, ident.1 = A, ident.2 = B,
                    min.pct = 0.1, logfc.threshold = 0.15, assay = "SCT")
degs$genes <- rownames(degs)
degs <- degs[order(degs$avg_log2FC, decreasing = T),]
wb <- createWorkbook()
addWorksheet(wb, sheet = paste0(A,".", B))
writeData(wb, sheet = paste0(A,".", B), degs)
saveWorkbook(wb, file=paste(rds@project.name, A, B, "DEGs.xlsx", sep = "-"), overwrite = TRUE)


## For all clust ####
A = "IL6"
B = "WT"
group = "geno"
Idents(rds) <- group
idents <- Idents(rds)
idents

rds <- PrepSCTFindMarkers(rds)
wb <- createWorkbook()
for(i in idents) {
  print(i)
  clusterA <- WhichCells(rds, idents = paste(A))
  clusterB <- WhichCells(rds, idents = paste(B))
  degs <- FindMarkers(rds, only.pos = F, ident.1 = clusterA, ident.2 = clusterB,
                      min.pct = 0.1, logfc.threshold = 0.15, assay = "SCT")
  degs$genes <- rownames(degs)
  degs <- degs[order(degs$avg_log2FC, decreasing = T),]
  addWorksheet(wb, sheet = i)
  writeData(wb, sheet = i, degs)
}
saveWorkbook(wb, file=paste(rds@project.name, A, B,  "DEGs.xlsx", sep = "-"), overwrite = TRUE)

# Subcluster ####
# For laptops/PCs override default variable size limit
n <- 50     # Let variables reach up to n GB
options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB

wkdir <- getwd()
original <- rds
Idents(original) <- original$basicAnno
idents <- levels(Idents(original))
idents
resRange <- seq(0.05, 0.2, by = 0.05)
alg <- 4
rnaPCs <- 30
adtPCs <- 20
algKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")

# Split cells from cluster back to separate objects by cartridge, rerun integration
for (cluster in idents) {
  set.seed(99)
  # Make new folder for each subcluster
  print(cluster)
  dir.create(paste0(wkdir, "/", cluster))
  setwd(paste0(wkdir,"/", cluster))
  rds <- subset(original, idents = cluster)
  list  <- SplitObject(rds, split.by = "orig.ident")

  ## RNA Normalization #####
  print("____________________________________________ RNA SCT nomalization start")
  list <- lapply(list, function(i) {
    DefaultAssay(i) <- 'RNA'
    i <- SCTransform(i, vars.to.regress = "percent.mt")
    return(i)
  })
  list
  print(paste(Sys.time(), "RNA SCT norm done"))
  gc()
  
  # RNA Integration ############## 
  print("____________________________________________ RNA integration start")
  features <- SelectIntegrationFeatures(list, nfeatures = 3000)
  list <- PrepSCTIntegration(list, anchor.features = features)
  list <- lapply(list, function(x) RunPCA(x, features = features, reduction.name = "rpca"))
  anchors <- FindIntegrationAnchors(list, anchor.features = features, dims = 1:rnaPCs,
                                    normalization.method = "SCT", reduction = "rpca")
  rds <- IntegrateData(anchors, normalization.method = "SCT", k.weight = 46)
  DefaultAssay(rds) <- "integrated"
  rds@project.name <- cluster
  # update rPCA reduction name
  rds[['subrpca']] <- rds[['pca']]
  gc()
  
  # Print PC Elbowplot
  pdf(paste0(rds@project.name,"-Elbow.pdf"))
  ElbowPlot(rds, reduction = "subrpca")
  dev.off()
  
  # Add old metadata
  rds[['rpca']] <- original[['rpca']]
  rds[['harmony.rna']] <- original[['harmony.rna']]
  rds[['rna.umap']] <- original[['rna.umap']]
  rds[['apca']] <- original[['apca']]
  rds[['adt.umap']] <- original[['adt.umap']]
  rds[['umap']] <- original[['umap']]
  rds[['wnn.umap']] <- original[['wnn.umap']]
  
  ## ADT CLR ####
  # print("____________________________________________ ADT CLR Normalization")
  # abseq <- rownames(rds@assays$ADT)
  # DefaultAssay(rds) <- 'ADT'
  # rds <- NormalizeData(rds, normalization.method = 'CLR')
  # rds <- ScaleData(rds)
  # DefaultAssay(rds) <- "SCT"
  
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
    positive_peak = list(ADT = "CD19"), 
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
  rds <- ScaleData(rds)
  rds <- RunPCA(rds, npcs = rnaPCs, reduction.name = "subrpca")
  VariableFeatures(rds@assays$ADT) <- rownames(rds@assays$ADT)
  rds <- RunPCA(rds, assay = "ADT", reduction.name = "subapca", npcs = adtPCs)
  
  ## Harmonize RNA ####
  # Batch effect correction
  rds <- RunHarmony(rds, group.by.vars = "orig.ident", 
                    reduction.use = "subrpca", reduction.save = "subharmony.rna")
  
  ### UMAP RNA ####
  # UMAP Without Harmony
  rds <- RunUMAP(rds, reduction = "subrpca", dims = 1:rnaPCs, reduction.name = "sub.rumap")
  # UMAP With Harmony
  rds <- RunUMAP(rds, reduction = "subharmony.rna", dims = 1:rnaPCs, reduction.name = "subharmony.rumap")
  
  # Print
  pdf(paste0(rds@project.name,"-PCA-Harmony.pdf"))
  Idents(rds) <- "orig.ident"
  print(DimPlot(rds, reduction = "sub.rumap"))
  print(DimPlot(rds, reduction = "subharmony.rumap"))
  print(DimPlot(rds, reduction = "subapaca"))
  print(DimPlot(rds, reduction = "subrpca"))
  dev.off()
  
  ## Multi-res Clust ####
  ### RNA ONLY ####
  print("_____________________ Starting RNA ONLY Clustering")
  set.seed(99)
  rds <- FindNeighbors(rds, dims = 1:rnaPCs, reduction = "subharmony.rna")
  
  # Multi-res Clustering
  pdf(paste0(rds@project.name, "-UMAP-RNA.pdf"))
  for(res in resRange) {
    print(res)
    rds <- FindClusters(rds, algorithm = alg, resolution = res)
    print(DimPlot(rds, reduction = "subharmony.rumap", label = TRUE, raster = F) +
            labs(title = paste0("RNA UMAP ", algKey[alg],": ", res)))
  }
  dev.off()
  
  print(paste0(Sys.time(), " -> RNA ONLY clustering done!"))
  
  ### RNA + ADT ####
  print("_____________________ Starting RNA & ADT Clustering")
  set.seed(99)
  rds <- FindMultiModalNeighbors(rds, reduction.list = list("subrpca", "subapca"),
                                 dims.list = list(1:rnaPCs, 1:adtPCs))
  rds <- RunUMAP(rds, nn.name= "weighted.nn", reduction.name = "subwnn.umap")
  # Multi-res Clustering
  pdf(paste0(rds@project.name, "-UMAP-RNAADT.pdf"))
  for(res in resRange) {
    print(res)
    rds <- FindClusters(rds, algorithm = alg, resolution = res)
    #, graph.name = "wsnn"
    print(DimPlot(rds, reduction= "subwnn.umap", label = T, raster = F) + 
            labs(title = paste0("ADT&RNA UMAP ", algKey[alg],": ", res)))
  }
  dev.off()
  
  ## SAVE!! ####
  saveRDS(rds, file= paste0(rds@project.name, "-RNAADT.rds"))
  
  print(paste0(Sys.time(), " -> RNA & ADT clustering done!"))
  
  ## QC FeatPlot ####
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "ADT_total")
  pdf(paste0(rds@project.name, "-QCUMAP.pdf"))
  print(FeaturePlot(rds, reduction = "subharmony.rumap", features = features, ncol = 2))
  print(FeaturePlot(rds, reduction = "subwnn.umap", features = features, ncol = 2))
  dev.off()
  
  ## Umap orig.ident ####
  pdf(paste(rds@project.name, "UMAP-batchEffectCheck.pdf", sep = "-"))
  print(DimPlot(rds, reduction = "subharmony.rumap", group.by = "orig.ident"))
  print(DimPlot(rds, reduction = "subharmony.rumap", group.by = "geno"))
  print(DimPlot(rds, reduction = "subwnn.umap", group.by = "orig.ident"))
  print(DimPlot(rds, reduction = "wsubnn.umap", group.by = "geno"))
  dev.off()
  
  ## Abseq Featplot #### 
  abseq <- rownames(rds@assays$ADT)
  
  # PRINT
  pdf(paste(rds@project.name, "AbsFeats-RNA.pdf", sep = "-"), width = length(abseq) / 5 * 4.8, height = length(abseq) / 5 * 5)
  rds@misc$umap <- "subharmony_rumap"
  DefaultAssay(rds) <-"ADT"
  plots <- FeaturePlot(rds, reduction = rds@misc$umap, features = abseq, ncol = 5) & 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(), 
          axis.ticks = element_blank())
  DefaultAssay(rds) <-"SCT"
  print(plots+ plot_annotation(title ="Abseq", theme = theme(plot.title = element_text(size =50))))
  dev.off()
  pdf(paste(rds@project.name, "AbsFeats-WNN.pdf", sep = "-"), width = length(abseq) / 5 * 4.8, height = length(abseq) / 5 * 5)
  rds@misc$umap <- "subwnn.umap"
  DefaultAssay(rds) <-"ADT"
  plots <- FeaturePlot(rds, reduction = rds@misc$umap, features = abseq, ncol = 5) & 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(), 
          axis.ticks = element_blank())
  DefaultAssay(rds) <-"SCT"
  print(plots+ plot_annotation(title ="Abseq", theme = theme(plot.title = element_text(size =50))))
  dev.off()
  
}

# Subcluster pt2 ####
## Pick res ####
rds$seurat_clusters<- rds$SCT_snn_res.0.05
rds@misc$umap <- "harmony_umap"
rds@project.name <- paste0("CMP-GMP Progenitor-RNA")
## VlnPlot ####
group <- "seurat_clusters"
Idents(rds) <- group
var2 <- "geno"

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

## Abseq dot ####
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
pdf(paste(rds@project.name, group, "AbseqDot.pdf", sep = "-"),width = 10, height = 6)
DotPlot(rds, features = abseq, col.min = 0, cols = "RdYlBu", group.by = group) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(hjust = 0)) + 
  ggtitle("AbSeq") + 
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5)
dev.off()
DefaultAssay(rds) <- "SCT"

## Abseq Ridge ####
DefaultAssay(rds) <- "ADT"
abseq <- rownames(rds@assays$ADT)
plots <- RidgePlot(rds, features = abseq, stack = T) + NoLegend()
pdf(paste0(rds@project.name, "-AbseqRidge.pdf"), width = 12, height = 5)
print(plots + labs(title = paste("Abseq Ridge Plot:", group)) + 
        theme(axis.text.y = element_text(hjust = 0))
)
dev.off()

DefaultAssay(rds) <- "SCT"

## Subset umap ####
split.by <- var2

pdf(paste(rds@project.name, "UMAPsubset.pdf", sep = "-"))
DimPlot(rds, reduction = rds@misc$umap, raster = F, group.by = group,split.by = split.by, ncol=2) + 
  NoLegend()
dev.off()

## Barplots ####
table <- data.frame(table(rds$seurat_clusters, rds[[var2]][,1]))
table$Percentage <- round((table$Freq / sum(table$Freq)) * 100, 1)

pdf(paste0(rds@project.name, "-Bar.pdf"), width = 5, height = 5)
ggplot(table, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_text(aes(label = Freq), # Show count % Pct
  #           vjust = -0.5, size =1) +
  labs(title = "Cell Type Distribution", x = "Cluster", y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
dev.off()

var1 <- group
data <- data.frame(table(rds[[var1]][,1], rds[[var2]][,1]))
colnames(data) <- c(var1, var2, "Freq")
data$Percentage <- round((data$Freq / sum(data$Freq)) * 100, 1)

# Write to csv file
write.csv(table(rds[[var1]][,1], rds[[var2]][,1]), paste0(rds@project.name,"-table.csv"))

# Stacked + percent
pdf(paste0(rds@project.name,"-BarPct.pdf"))
print(ggplot(data, aes(fill=geno, y=Percentage, x=seurat_clusters)) +   #EDIT
        geom_bar(position="fill", stat="identity") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Cell Composition"))
dev.off()


## DEG Subclusters####
group <- "condition"
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

### Genos ####
degs <- NULL
wb <-NULL
group <- "geno"
Idents(rds) <- group
comparisons <- list(list("IL6","WT"), list("GCSF", "WT"), list("DKO", "WT"))
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

### Filter DEGs ####
# by LogFC, gate out upregs if pct1 < 0.5
### IF BLANK SHEETS MANUALLY DELETE!!!!!! ####
input_file_list <- c(paste(rds@project.name, A, B, "DEGs", sep = "-"))

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



### Dotplot top 10 DEGs ####
# FILTER OUT SMALL CLUSTERS
# n <- 50
# rds$clust.cond <- paste(rds$seurat_clusters, rds$geno)
# Idents(rds) <- rds$clust.cond
# clustSize <- table(rds$clust.cond)
# keep  <- names(clustSize[clustSize >= n])
# sub <-subset(rds, idents = keep)
#### by cluster ####
top <- 10
input_file <- paste0("EndoFocus-anno.geno-DEGs")
sheet_names <- excel_sheets(paste0(input_file, "-filtered.xlsx"))
genes <- c()
for (sheet in sheet_names) {
  print(sheet)
  # Read the sheet
  degs <- read.xlsx(paste0(input_file, "-filtered.xlsx"), sheet = sheet)
  genes <- append(genes, degs$gene[1:top])
}
genes <- unique(genes[!is.na(genes)])
genes <- genes[!grepl("^ENS", genes)]
genes <- genes[!grepl("^LINC", genes)]
pdf(paste0(rds@project.name, "-seurat_clusters-DotTop", top,".pdf"), width =22, height =4)
print(DotPlot(rds, features = unique(genes), cols = "RdYlBu", col.min = 0, dot.scale = 5, group.by = "anno.geno") +
        ggtitle(paste(rds@project.name, "Top 10 per cluster")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) +
        geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
        geom_hline(yintercept = seq(n+0.5, length(unique(rds$anno.geno)) - 0.5, by = n), color = "black")
)
dev.off()

#### by condition ####
# FILTER OUT SMALL CLUSTERS
# n <- 50
# rds$cond.clust <- paste(rds$geno, rds$seurat_clusters)
# Idents(rds) <- "cond.clust"
# clustSize <- table(rds$cond.clust)
# keep  <- names(clustSize[clustSize >= n])
# sub <-subset(rds, idents = keep)

input_file <-paste0(rds@project.name, "-geno-DEGs")
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
        ggtitle(paste(rds@project.name, "Top 10 per genotype")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) +
        geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5))
dev.off()

# Pathway Geno ####
# List of Enrichr databases: https://maayanlab.cloud/Enrichr/#libraries
databaseList <- c("KEGG_2019_Mouse", 
                  "GO_Biological_Process_2021",
                  "Mouse_Gene_Atlas", 
                  "MSigDB_Hallmark_2020",
                  "Reactome_Pathways_2024",
                  "WikiPathways_2024_Mouse")

Idents(rds) <- group
comparisons <- list(list("IL6","WT"), list("GCSF", "WT"), list("DKO", "WT"))
for (pair in comparisons){
  A = pair[[1]]
  B = pair[[2]]
  pdf(paste0(rds@project.name, "-Pathway-", A, "v", B, ".pdf"), width = 10, height = 6)
  for(db in databaseList){
    print(DEenrichRPlot(rds, ident.1 = A, ident.2 = B, max.genes = 100,
                        logfc.threshold = 0.25, p.val.cutoff = 0.05, num.pathway = 10,
                        enrich.database = db))
  }
  dev.off() 
}


# Subset rds by features ####
rds <- micro
keep <- c("CD45", "CD11b", "CD115", "F4.80", "CD11c")
rds[['ADT']] <- subset(micro[['ADT']], features= keep)