set.seed(99)
library(openxlsx)
library(Seurat)
library(RColorBrewer)
library(edgeR)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(scDblFinder)
library(dplyr)
library(SingleCellExperiment)
library(harmony)
library(sctransform)
library(purrr)
library(leidenAlg)
library(Seurat)
setwd("/mnt/BioAdHoc/Groups/Collaborators/ben.croker/nomid")

test <- "bimod"
A = "Brain_NM"
B = "Brain_WT"
group = "tissue.strain"
Idents(rds) <- rds$annotations
idents <- levels(Idents(rds))
idents

Idents(rds) <- rds$tissue.strain

# DEG for whole nomid unST
wb <- createWorkbook()
for(i in idents) {
  print(i)
  clusterA <- WhichCells(rds, idents = paste(i, A))
  clusterB <- WhichCells(rds, idents = paste(i, B))
  degs <- FindMarkers(rds, test.use = test, only.pos = F, ident.1 = clusterA, ident.2 = clusterB,
                      min.pct = 0.1, logfc.threshold = 0.15, assay = "SCT")
  degs$genes <- rownames(degs)
  degs <- degs[order(degs$avg_log2FC, decreasing = T),]
  addWorksheet(wb, sheet = paste0(A,".", B))
  writeData(wb, sheet = paste0(A,".", B), degs)
}
saveWorkbook(wb, file=paste(rds@project.name, A, B, test, "DEGs.xlsx", sep = "-"), overwrite = TRUE)


# Custom Dotplot lines ####
custom <- c(4.5 , 8.5, 12.5, 16.5, 20.5, 24.5, 26.5, 29.5, 33.5, 37.5 ,40.5, 44.5 ,47.5, 49.5, 52.5, 54.5, 58.5, 60.5,
            64.5, 68.5, 72.5, 76.5, 80.5, 83.5, 87.5, 91.5, 95.5)
listName <- "AllGOI"
group <- "anno.tissue.strain"
n = 4
Idents(rds) <- rds[[group]][,1]
DefaultAssay(rds) <- "SCT"
genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)

# PRINT
pdf(paste(rds@project.name, listName, group, "Dot.pdf", sep = "-"), height =10, width = 20)
for (col in names){
  print(col)
  print(DotPlot(rds, features = genes[col], cols = "RdYlBu", col.min = 0, dot.scale = 5) +
          coord_flip() + 
          ggtitle(col) +
          #geom_hline(yintercept = seq(n+0.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black") +
          geom_hline(yintercept = custom, color = "black") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)))
}
dev.off()

# Erythroid custom lines ####
custom <- c(4.5,  8.5, 12.5, 16.5, 20.5, 24.5, 28.5, 32.5, 35.5, 39.5, 43.5, 47.5)
listName <- "dot2"
group <- "cluster.tissue.strain"
Idents(rds) <- rds[[group]][,1]
DefaultAssay(rds) <- "SCT"
genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)

# PRINT
pdf(paste(rds@project.name, listName, group, "Dot.pdf", sep = "-"), height =10, width = 60)
for (col in names){
  print(col)
  print(DotPlot(rds, features = genes[col], cols = "RdYlBu", col.min = 0, dot.scale = 5) +
          #coord_flip() + 
          ggtitle(col) +
          geom_hline(yintercept = custom, color = "black") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)))
}
dev.off()

# subcluster to whole
Idents(rds) <- rds$annotations
Idents(neuts) <- neuts$lineage
table(Idents(rds))
table(Idents(neuts))

sample <- WhichCells(rds, idents = "Immature Neutrophil")

Idents(rds, cells = sample) <- "Immature Neut 2"

DimPlot(rds, reduction = "rna.umap", label = T) + NoLegend()

rds$annotations <- Idents(rds)

# Subclust to Whole ####
rds <- nomid
Idents(rds) <- rds$basicAnno
Idents(micro) <- micro$annotations
Idents(neuts) <- neuts$annotations
Idents(endo) <- endo$subclustwsnn_res.0.2

# Micro
idents <- levels(Idents(micro))
idents
for(i in idents){
  select.cells <- WhichCells(micro, idents = i)
  Idents(rds, cells = select.cells) <- i
}

# Endo
idents <- levels(Idents(endo))
idents
for(i in idents){
  select.cells <- WhichCells(endo, idents = i)
  Idents(rds, cells = select.cells) <- paste("Endothelial", i)
}

# Neuts
idents <- levels(Idents(neuts))
idents
for(i in idents){
  select.cells <- WhichCells(neuts, idents = i)
  Idents(rds, cells = select.cells) <- paste(i)
}

# alphabetize the cell types
Idents(rds) <- factor(Idents(object = rds), levels = sort(levels(rds)))
rds$subcluster <- Idents(rds)

#PRINT
pdf(paste0(rds@project.name, "-subclustered.pdf"), width = 12, height = 8)
DimPlot(rds, reduction = rds@misc$umap, raster = F, label = T, label.size = 3)
dev.off()


# DEGs w/subclusters ####
# Remove residuals to force recalculation
DefaultAssay(rds) <- "integrated"
# Set the default assay to SCT
DefaultAssay(rds) <- "SCT"

# Access the UMI counts matrix
umi_counts <- GetAssayData(rds, slot = "counts")

# Calculate the minimum UMI per gene across all cells
min_umi_per_gene <- rowSums(umi_counts > 0)

# Get the minimum value
min_umi <- min(min_umi_per_gene)
min_umi

rds@assays$SCT@misc$residuals

rds@assays$integrated@misc$min_umi <- 1
PrepSCTFindMarkers(rds)

# Split rds for upload ####
types <- rownames(table(Idents(rds)))
types

for (cluster in types){
  sub <- subset(rds, idents = cluster)
  saveRDS(sub, paste0("Deaddying.rds"))
}

# Abseq focus Dot CLR ####
pdf(paste0(rds@project.name, "-AbseqDot.pdf"), width = 12, height = 13)
DotPlot(rds, features = nomidAbseq, col.min = 0, cols = "RdYlBu") + 
  ggtitle(paste("AbSeq CLR", group)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) + 
  geom_hline(yintercept = seq(n-1.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(n+0.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black")
dev.off()

# Abseq subset ####
focus<- c("Endothelial Cell Spleen WT", "Endothelial Cell Spleen NM", "Endothelial Cell Brain NM", "Endothelial Cell Brain WT",
          "M1 Microglia Brain NM", "M1 Microglia Brain WT", "M1 Microglia Spleen NM",
          "M1 Microglia Spleen WT", "M2 Microglia Brain NM", "M2 Microglia Brain WT", 
          "M2 Microglia Spleen NM", "M2 Microglia Spleen WT", "Mature Neutrophil Brain NM", 
          "Mature Neutrophil Brain WT", "Mature Neutrophil Spleen NM",  "Mature Neutrophil Spleen WT",
          "NeP Brain NM", "NeP Brain WT", "NeP Spleen NM", "NeP Spleen WT")
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
pdf(paste(rds@project.name, "AbseqDot-subset.pdf", sep = "-"),width = 10, height =7)
DotPlot(rds, features = abseq, col.min = 0, cols = "RdYlBu", idents = focus) + 
  #geom_hline(yintercept = seq(n+0.5, length(focus) - 0.5, by = n), color = "black") +
  geom_hline(yintercept = seq(n-1.5, length(focus) - 0.5, by = n), color = "black", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("AbSeq")
dev.off()
DefaultAssay(rds) <- "SCT"

n <- 2
focus<- c("Endothelial Cell Spleen WT", "Endothelial Cell Spleen NM", "Endothelial Cell Brain NM", "Endothelial Cell Brain WT")
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
pdf(paste("Endothelial-AbseqDot-subset.pdf", sep = "-"),width = 10, height =3.5)
DotPlot(rds, features = abseq, col.min = 0, cols = "RdYlBu", idents = focus) + 
  geom_hline(yintercept = seq(n+0.5, length(focus) - 0.5, by = n), color = "black") +
  geom_hline(yintercept = seq(n-0.5, length(focus) - 0.5, by = n), color = "black", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Endothelial Cell AbSeq")
dev.off()
DefaultAssay(rds) <- "SCT"

n <- 4
focus<- c("M1 Microglia Brain NM", "M1 Microglia Brain WT", "M1 Microglia Spleen NM",
          "M1 Microglia Spleen WT", "M2 Microglia Brain NM", "M2 Microglia Brain WT", 
          "M2 Microglia Spleen NM", "M2 Microglia Spleen WT")
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
pdf(paste("Microglial-AbseqDot-subset.pdf", sep = "-"),width = 10, height =4.5)
DotPlot(rds, features = abseq, col.min = 0, cols = "RdYlBu", idents = focus) + 
  geom_hline(yintercept = seq(n+0.5, length(focus) - 0.5, by = n), color = "black") +
  geom_hline(yintercept = seq(n-1.5, length(focus) - 0.5, by = n), color = "black", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Microglia AbSeq")
dev.off()
DefaultAssay(rds) <- "SCT"

focus<- c("Mature Neutrophil Brain NM", "Mature Neutrophil Brain WT", "Mature Neutrophil Spleen NM",  "Mature Neutrophil Spleen WT",
          "NeP Brain NM", "NeP Brain WT", "NeP Spleen NM", "NeP Spleen WT")
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
pdf(paste("Neutrophils-AbseqDot-subset.pdf", sep = "-"),width = 10, height =4.5)
DotPlot(rds, features = abseq, col.min = 0, cols = "RdYlBu", idents = focus) + 
  geom_hline(yintercept = seq(n+0.5, length(focus) - 0.5, by = n), color = "black") +
  geom_hline(yintercept = seq(n-1.5, length(focus) - 0.5, by = n), color = "black", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Neutrophils AbSeq")
dev.off()
DefaultAssay(rds) <- "SCT"

# Abseq with subsets ####
# Endo
idents <- levels(Idents(rds))
focus<- grep(pattern = "Endothelial", x = idents, value = T)
focus<- focus[-length(focus)]
focus
abseq <- rownames(rds[["ADT"]])
#abseq <- c("CD105",  "CD71", "CD184")
DefaultAssay(rds) <- "ADT"
pdf(paste("Endothelial-AbseqDot-subclust.pdf", sep = "-"),width = 4.5, height =4.5)
DotPlot(rds, features = abseq, col.min = 0, cols = "RdYlBu", idents = focus) + 
  #geom_hline(yintercept = seq(n+0.5, length(focus) - 0.5, by = n), color = "black") +
  #geom_hline(yintercept = seq(n-1.5, length(focus) - 0.5, by = n), color = "black", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Endothelial Cells AbSeq")
dev.off()
DefaultAssay(rds) <- "SCT"


# Neut
idents <- levels(Idents(rds))
focus<- grep(pattern = "Neut", x = idents, value = T)
focus <- c(focus, "NeP", "Progenitor", "Eosinophil")
abseq <- rownames(rds[["ADT"]])
abseq <- c("CD11b",  "CD16/32", "CD45", "CD62L", "CXCR2", "Ly6G", "SiglecF")
DefaultAssay(rds) <- "ADT"
pdf(paste("Neutrophils-AbseqDot-subclust.pdf", sep = "-"),width = 6, height =4.5)
DotPlot(rds, features = abseq, col.min = 0, cols = "RdYlBu", idents = focus) + 
  #geom_hline(yintercept = seq(n+0.5, length(focus) - 0.5, by = n), color = "black") +
  #geom_hline(yintercept = seq(n-1.5, length(focus) - 0.5, by = n), color = "black", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Neutrophils AbSeq")
dev.off()
DefaultAssay(rds) <- "SCT"

# DEG ####
group <- "Microglia"
focus<- c("M1 Microglia", "M1 Microglia", "M1 Microglia",
          "M1 Microglia", "M2 Microglia", "M2 Microglia", 
          "M2 Microglia", "M2 Microglia")
focus
idents <- levels(Idents(rds))
rds <- PrepSCTFindMarkers(rds)

## DEG.xlsx -> Volcano ####
#input_file_list <- c("T Cell-GCSF-IL6-DEGs-wilcox", "Monocyte-GCSF-IL6-DEGs-wilcox")

for (input_file in input_file_list){
  print(input_file)
  sheet_names <- excel_sheets(paste0(input_file, ".xlsx"))
  pdf(paste0(input_file, "-Volcano.pdf"))
  
  for (sheet in sheet_names){
    print(sheet)
    degs <- read.xlsx(paste0(input_file, ".xlsx"), sheet = sheet)
    print(EnhancedVolcano(degs, lab = degs$gene, x = "avg_log2FC", y = "p_val_adj",
                          pCutoff = 0.05, FCcutoff = 0.5,
                          title = sheet, subtitle = ""))
  }
  dev.off()
}

# Pseudobulk ####
# Add Metadata
rds <- AddMetaData(rds, metadata = rds$annotations, col.name = "annotations")
head(rds@meta.data$annotations)
table(rds@meta.data$annotations)

rds <- AddMetaData(rds, metadata = rds$geno, col.name = "geno")
head(rds@meta.data$geno)
table(rds@meta.data$geno)

# Neuts
Idents(rds) <- rds$annotations
pseudo  <- AggregateExpression(rds, assays= "SCT", return.seurat = T, 
                               group.by = c("annotations", "geno"))
table(pseudo$orig.ident)
# Metadata is lost when AggregateExpression is run
# Add back metadata
meta <- data.frame(annotations = sapply(strsplit(Cells(pseudo), "_"), `[`, 1),
                   geno = sapply(strsplit(Cells(pseudo), "_"), `[`,2),
                   row.names = Cells(pseudo)
)
pseudo <- AddMetaData(pseudo, meta)

# Make new identifier col
pseudo$anno.geno <- paste(pseudo$annotations,pseudo$geno)
Idents(pseudo) <- pseudo$anno.geno
table(Idents(pseudo))

pseudoDeg <- FindAllMarkers(pseudo, test.use = "DESeq2", group.by = "")
pseudoDeg <- FindMarkers(pseudo, test.use = "DESeq2",
                         ident.1 = "Brain Mature Neutrophil",
                         ident.2 = "Spleen Mature Neutrophil",)

head(pseudoDeg, n = 15)
## Just WT vs NM ####
expr <- GetAssayData(pseudo, slot = "data")
ident1 <- expr[, grep("WT$", colnames(expr))]
ident2 <- expr[, grep("NM$", colnames(expr))]

# Compute log2 fold change
logFC <- rowMeans(ident1) -  rowMeans(ident2)
logFC <- sort(logFC, decreasing = T)
head(logFC)

# Pathway Enrichment ####
# clusterProfiler
install.packages("clusterProfiler")
library(org.Mm.eg.db)

input_file <- "Nomid-WNN-Neutrophils-DEGs-filtered"
print(input_file)
sheet_names <- excel_sheets(paste0(input_file, ".xlsx"))

for (sheet in sheet_names){
  print(sheet)
  degs <- read.xlsx(paste0(input_file, ".xlsx"), sheet = sheet)
}

# VlnPlot tiss.geno.basic####
pdf(paste0(rds@project.name, "-QCVln-", group, ".pdf"), width = 15, height = 6)
for (feature in features){
  print(VlnPlot(rds, features = feature, pt.size = 0, group.by = group) + 
          geom_vline(xintercept = seq(n+0.5, length(unique(Idents(rds))) - 0.5, by = n), linetype = "dashed", color = "black") + 
          NoLegend())
}
dev.off()

# Abseq Edited Ridge ####
group <-"tiss.geno.sub"
Idents(rds) <- group
DefaultAssay(rds) <- "ADT"
plots <- RidgePlot(rds, features = abseq, stack = T) + NoLegend()
# PRINT
pdf(paste0(rds@project.name, "-AbseqRidge-", group, "edit.pdf"), width = 20, height = 40)
print(plots + labs(title = paste("Abseq Ridge Plot:", group)) + 
        theme(axis.text.y = element_text(hjust = 0))
)
dev.off()
DefaultAssay(rds) <- "SCT"

#_______________________
group  <- "basicAnno"
Idents(rds) <- group
DefaultAssay(rds) <- "ADT"
abseq <- c("CD105", "CD117", "CD11a", "CD11b", "CD162", "CD16/32",
           "CD184", "CD19", "CD45", "CD48", "CD62L", "CD71", "CXCR2",
           "CD150", "F4/80", "Ly6A/E", "Ly6G", "SiglecF", "TCRb")
plots <- RidgePlot(rds, features = abseq, stack = T) + NoLegend()
# PRINT
pdf(paste0(rds@project.name, "-AbseqRidge-", group, "-edit.pdf"), width = 20, height = 10)
print(plots + labs(title = paste("Abseq Ridge Plot:", group)) + 
        theme(axis.text.y = element_text(hjust = 0))
)
dev.off()

#____________________________
group  <- "tiss.geno.basic"
Idents(rds) <- group
DefaultAssay(rds) <- "ADT"
plots <- RidgePlot(rds, features = abseq, stack = T) + NoLegend()
# PRINT
pdf(paste0(rds@project.name, "-AbseqRidge-", group, "-edit.pdf"), width = 20, height = 10)
print(plots + labs(title = paste("Abseq Ridge Plot:", group)) + 
        theme(axis.text.y = element_text(hjust = 0))
)
dev.off()
DefaultAssay(rds) <- "SCT"

#_____________________________
group  <- "subcluster"
Idents(rds) <- group
DefaultAssay(rds) <- "ADT"
plots <- RidgePlot(rds, features = abseq, stack = T) + NoLegend()
# PRINT
pdf(paste0(rds@project.name, "-AbseqRidge-", group, "-edit.pdf"), width = 20, height = 30)
print(plots + labs(title = paste("Abseq Ridge Plot:", group)) + 
        theme(axis.text.y = element_text(hjust = 0))
)
dev.off()
DefaultAssay(rds) <- "SCT"

# ADTnorm ####
BiocManager::install("flowCore")
remotes::install_github("yezhengSTAT/ADTnorm", build_vignettes = F)
library(ADTnorm)

# Fix no "/" in abseq names
newAbseq <- c("CD105" , "CD115" , "CD117",  "CD11a",   "CD11b" ,  "CD11c" ,  "CD162" ,  "CD16.32",
              "CD184",   "CD19",    "CD335",   "CD41",    "CD45",    "CD48" ,   "CD62L",   "CD71",
              "CXCR2",   "Clec7a",  "CD150" ,  "F4.80",   "Ly6A.E",  "Ly6G" ,   "NK1.1" ,  "SiglecF",
              "TCRb",    "TER119")

temp <- rds@assays$ADT@counts
rownames(temp) <- newAbseq
rds[["ADT"]] <- CreateAssayObject(counts = temp)

## nomidAbseq ####
# Select markers useful for nomid
nomidAbseq <- c("CD105", "CD117", "CD11a", "CD11b", "CD162", "CD16.32",
                "CD184", "CD19", "CD45", "CD48", "CD62L", "CD71", "CXCR2",
                "CD150", "F4.80", "Ly6A.E", "Ly6G", "SiglecF", "TCRb")


group <-"basicAnno"
#group <- "tiss.geno.basic"
Idents(rds)<- group

# Get raw ADT counts with cell ID rows, Abseq features cols
cell_x_adt <- t(rds@assays$ADT@counts)
head(cell_x_adt)

# Get desired idents 
cell_x_feature <- data.frame(rds@active.ident)

# Make sure there is a 'sample' column set to an ident
colnames(cell_x_feature) <- 'sample'
head(cell_x_feature)

cell_x_adt_norm = ADTnorm(
  cell_x_adt = cell_x_adt, 
  cell_x_feature = cell_x_feature, 
  save_outpath = getwd(), 
  study_name = rds@project.name,
  marker_to_process = newAbseq,
  bimodal_marker = NULL,             # default NULL: try different settings to find bimodal peaks
  trimodal_marker = c("CD45"),       # CD4 and CD45RA tend to have 3 peaks
  # setting the CD3 uni-peak of buus_2021_T study to positive peak if only one peak is detected for CD3 marker
  # positive_peak = list(ADT = "CD3", sample = "buus_2021_T"), 
  positive_peak = list(ADT = "CD19", sample = "B Cell"), 
  brewer_palettes = "Dark2",
  save_fig = TRUE,
  target_landmark_location = "fixed",
  shoulder_valley = T,               # Look for "shoulder" as pos peak (technical variation -> no clear separation b/w neg/pos)
  #multi_sample_per_batch = T,        # Omit aligning the one pos peak
  #customize_landmark = T             # Manual adjustment UI 
)

## Put ADTNorm matrix back into rds ####
rds@assays$ADT@data <- t(cell_x_adt_norm)

## VlnPlot ####
pdf(paste0(rds@project.name,"-adtNorm-VlnPlot.pdf"))
VlnPlot(rds, features = nomidAbseq, assay = "ADT")
dev.off()

## RidgePlot ####
group <- "annotations"
Idents(rds) <- group
DefaultAssay(rds) <- "ADT"
plots <- RidgePlot(rds, features = nomidAbseq, stack = T) + NoLegend()
# PRINT
pdf(paste0(rds@project.name, "-AbseqRidge-", group, "-ADTNorm.pdf"), width = 20, height = 10)
print(plots + labs(title = paste("ADTNorm Ridge Plot:", group)) + 
        theme(axis.text.y = element_text(hjust = 0))
)
dev.off()

## FeatPlot ####
DefaultAssay(rds) <-"ADT"
plots <- FeaturePlot(rds, reduction = rds@misc$umap, features = nomidAbseq, ncol = 5) & 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank())
DefaultAssay(rds) <-"SCT"

# PRINT
pdf(paste(rds@project.name, "AbsFeat-UMAP.pdf", sep = "-"), width = length(nomidAbseq) / 5 * 4.8, height = length(nomidAbseq) / 5 * 5)
print(plots+ plot_annotation(title ="Abseq", theme = theme(plot.title = element_text(size =50))))
dev.off()


# Dotplot
rds$tiss.geno.basic <- paste(rds$basicAnno, rds$tiss.geno)
group <- "tiss.geno.basic"
Idents(rds)<- group
# Alphabetize
Idents(rds) <- factor(Idents(object = rds), levels = sort(levels(rds)))
rds$tiss.geno.basic <- Idents(rds)

n = 4
pdf(paste0(rds@project.name, "-AbseqDot-", group, "-edit-ADTNorm.pdf"), width = 12, height = 13)
DotPlot(rds, features = nomidAbseq, col.min = 0, cols = "RdYlBu") + 
  ggtitle(paste("AbSeq ADTNorm", group)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) + 
  geom_hline(yintercept = seq(n-1.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black", linetype = "dashed") +
  geom_hline(yintercept = seq(n+0.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black")
dev.off()


## Rerun WNN ####
alg <- 4
resRange <- seq(0.2, 0.4, by = 0.1)
rnaPCs <- 30
adtPCs <- 26
algKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")

DefaultAssay(rds)<-"ADT"
rds <- FindVariableFeatures(rds, features = rownames(rds[["ADT"]]))
rds <- ScaleData(rds, verbose = FALSE)
rds <- RunPCA(rds, npcs = length(rownames(rds@assays$ADT)), reduction.name = 'apca', 
              verbose = FALSE, approx = FALSE)
rds <- FindMultiModalNeighbors(rds, reduction.list = list("harmony.rna", "apca"),
                               dims.list = list(1:rnaPCs, 1:adtPCs), modality.weight.name = "integrated.weight")
rds <- RunUMAP(rds, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# Save metadata 
rds@project.name <- paste0(rds@project.name, "-adtNorm")
rds@misc$clustAlg <- algKey[alg]
rds@misc$umap <- "wnn.umap"

# DoubletFinder
# library(scDbFinder)
# sce <- as.SingleCellExperiment(rds, assay = "SCT")
# sce <- scDblFinder(sce, samples = "orig.ident")
# 
# rds$scDblFinder.class <- sce$scDblFinder.class
# rds$scDblFinder.class <- sce$scDblFinder.score
# 
# rds <- AddMetaData(rds, metadata=sce$scDblFinder.score, col.name='scDblFinder_score')
# rds <- AddMetaData(rds, metadata=sce$scDblFinder.class, col.name='scDblFinder_class')
# 
# pdf(paste0(project, "-RNA_dbFinderUMAP.pdf"))
# DimPlot(rds, reduction="rna.umap", raster = F, group.by="scDblFinder_class", cols=c('grey', 'red'), order=TRUE)
# FeaturePlot(rds, reduction="rna.umap", features='scDblFinder_score', raster = F, order=TRUE)
# dev.off()
# 
# Save summary
# write.csv(data.frame(prop.table(table(rds$scDblFinder_class))*100), paste0(project, "-DblFinder.csv"), row.names = F)
# write.table(data.frame(Class = "Doublet", 
#                        Statistic = names(summary(rds$nCount_RNA[rds$scDblFinder_class=='doublet'])),
#                        Value = as.numeric(summary(rds$nCount_RNA[rds$scDblFinder_class=='doublet']))),
#             paste0(project, "-DblFinder.csv"), sep = ",", row.names = F, col.names = F, append = T)
# write.table(data.frame(Class = "Singlet", 
#                        Statistic = names(summary(rds$nCount_RNA[rds$scDblFinder_class=='singlet'])),
#                        Value = as.numeric(summary(rds$nCount_RNA[rds$scDblFinder_class=='singlet']))),
#             paste0(project, "-DblFinder.csv"), sep = ",", row.names = F, col.names = F, append = T)
# # Save
# saveRDS(rds, paste0(project, "-DbFinder.rds"))
# Idents(rds) <- rds$scDblFinder_class
# rds <- subset(rds, idents= "singlet")
# 
# Multi-res Clustering
pdf(paste0(rds@project.name, "-WNNUMAP.pdf"))
for(res in resRange) {
  print(res)
  rds <- FindClusters(rds, algorithm = alg, resolution = res, graph.name = "wsnn")
  print(DimPlot(rds, reduction= "wnn.umap", label = T, raster = F) + 
          labs(title = paste0("ADT(adtNorm) + RNA UMAP ", algKey[alg],": ", res)))
}
dev.off()
saveRDS(rds, file= paste0(rds@project.name, "-WNN.rds"))

# ADTnorm UMAPS ####
pdf(paste0(rds@project.name, "-WNNUMAP-ANNO.pdf"), width = 12, height = 10)
DimPlot(rds, reduction = "wnn.umap",group.by = "annotations", raster = F, label = T) +
  # Make legend 1 column and adjust dot size
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))
dev.off()


# Custom order dotplots ####
## Neuts ####
custom <- c("NeP",
            "Pre-Neutrophil",
            "Immature Neutrophil",
            "Brain Mature Neutrophil",
            "Spleen Mature Neutrophil",
            "Aged Neutrophil",
            "Eosinophil"
            )
# order the cell types
Idents(rds) <- factor(Idents(rds), levels = rev(custom))
rds$annotations <- Idents(rds)

rds$anno.geno <-  paste(rds$annotations, rds$geno)
Idents(rds)<- "anno.geno"
levels(Idents(rds))
custom <- c("NeP NM",
            "NeP WT",
            "Pre-Neutrophil NM",
            "Pre-Neutrophil WT",           
            "Immature Neutrophil NM",
            "Immature Neutrophil WT",
            "Spleen Mature Neutrophil NM",
            "Spleen Mature Neutrophil WT", 
            "Brain Mature Neutrophil NM",
            "Brain Mature Neutrophil WT",
            "Aged Neutrophil NM",
            "Aged Neutrophil WT",
            "Eosinophil NM",
            "Eosinophil WT"
            )
# order the cell types
Idents(rds) <- factor(Idents(rds), levels = rev(custom))
rds$anno.geno <- Idents(rds)

## Endo ####
rds$anno.geno <- paste(rds$annotations, rds$geno)
Idents(rds) <- "anno.geno"
levels(Idents(rds))
custom <- c("Endothelial Healthy NM",
            "Endothelial Healthy WT",   
            "Endothelial Inflammed NM", 
            "Endothelial Inflammed WT", 
            "Endothelial Mitotic NM",
            "Endothelial Mitotic WT",
            "Endothelial Artery NM",    
            "Endothelial Artery WT"
            )
# order the cell types
Idents(rds) <- factor(Idents(rds), levels = rev(custom))
rds$anno.geno <- Idents(rds)

## Micro ####
custom <- c("Microglial Progenitor",
            "M1 Microglia",
            "M2 Microglia",
            "Microglia Endothelial",
            "Vascular Endothelial"
            )
# order the cell types
Idents(rds) <- factor(Idents(rds), levels = rev(custom))
rds$annotations <- Idents(rds)


rds$anno.geno <- paste(rds$annotations, rds$geno)
Idents(rds) <- "anno.geno"
levels(Idents(rds))
custom <- c("Microglial Progenitor NM",
            "Microglial Progenitor WT",
            "M1 Microglia NM",
            "M1 Microglia WT",
            "M2 Microglia NM",
            "M2 Microglia WT",
            "Microglia Endothelial NM",
            "Microglia Endothelial WT",
            "Vascular Endothelial NM",
            "Vascular Endothelial WT"
            )
# order the cell types
Idents(rds) <- factor(Idents(rds), levels = rev(custom))
rds$anno.geno <- Idents(rds)


# CellChat ####
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
## Compare Cellchat ####
c1 <- readRDS("Nomid-WNN-adtNorm-CC-Brain NM.rds")
c2 <- readRDS("Nomid-WNN-adtNorm-CC-Brain WT.rds")
c3 <- readRDS("Nomid-WNN-adtNorm-CC-Spleen NM.rds")
c4 <- readRDS("Nomid-WNN-adtNorm-CC-Spleen WT.rds")

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