# install.packages('devtools')
# BiocManager::install('multtest')
# install.packages('metap')
# devtools::install_github('immunogenomics/presto')
set.seed(99)
library(dplyr)
if(!require("scCustomize", quietly = T)) install.packages("scCustomize")
if(!require("ComplexHeatmap", quietly = T)) BiocManager::install("ComplexHeatmap")
n
if(!require("openxlsx", quietly = T)) install.packages("openxlsx")
if(!require("readxl", quietly = T)) install.packages("readxl")
if(!require("EnhancedVolcano", quietly = T)) BiocManager::install("EnhancedVolcano")
n
if (!require("harmony", quietly = TRUE)) BiocManager::install("harmony")
n
install.packages("devtools")
devtools::install_github('immunogenomics/presto')
library(openxlsx)
library(readxl)
library(patchwork)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(EnhancedVolcano)
display.brewer.all(colorblindFriendly = TRUE)

# Compression of rds ####
library(qs)

# Save
qsave(rds, paste0(rds@project.name, ".qs"))

# Load
rds <- qread(paste0(rds@project.name, ".qs"))

#  QC post Cluster #####
# Cells/cluster table
write.csv(table(Idents(rds)), file=paste(rds@project.name, "UMAP.csv", sep = "-"))

features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "ADT_total")
rds$ADT_total <- colSums(rds@assays$ADT@data)

## FeatPlot ####
pdf(paste0(rds@project.name, "-QCUMAP.pdf"))
FeaturePlot(rds, reduction = rds@misc$umap, features = features, ncol = 2)
dev.off()

## VlnPlot ####
group <- "annotations"
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

# UMAP subset #####
## WhichCells way #####
group <- "geno"

Idents(rds) <- rds[[group]][,1]
idents <- levels(Idents(rds))
plots <- lapply(idents, function(i){
  print(i)
  cells <- WhichCells(rds, idents = i)
  DimPlot(rds, reduction =  rds@misc$umap, cells.highlight = cells, raster=F) +
    labs(title = paste("UMAP",i)) + NoLegend()
})

# PRINT
pdf(paste(rds@project.name, "UMAPsubset-WhichCells.pdf", sep = "-"),width = length(idents) / 3 * 6, height = length(idents) / 3.5 * 8)
print(wrap_plots(plots, ncol = 2))
dev.off()

## Subsets way #####
group <- "annotations"
split.by <- "condition"
Idents(rds) <- group

pdf(paste(rds@project.name, "UMAPsubset.pdf", sep = "-"), width = 8, height =5)
DimPlot(rds, reduction = "harmony_umap", raster = F, group.by = group,split.by = split.by, ncol=2) + 
  NoLegend()
dev.off()

# ABSEQ ################################# 
## DotPlot ####
group <-"annotations"
Idents(rds) <- group
n = 2
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
pdf(paste(rds@project.name, group, "AbseqDot.pdf", sep = "-"),width = 10, height = 6)
DotPlot(rds, features = abseq, col.min = 0, cols = "RdYlBu", group.by = group) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(hjust = 0)) + 
  ggtitle("AbSeq") +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5)
# coord_flip()
# geom_hline(yintercept = seq(n+0.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black")
#geom_hline(yintercept =seq(n-floor(n/2)+ 0.5, length(unique(Idents(rds))) - 0.5, by = n), linetype = "dashed")
#geom_hline(yintercept = custom, color = "black")
dev.off()
DefaultAssay(rds) <- "SCT"

#### scCustom ####
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
colors <- colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))
pdf(paste(rds@project.name, "AbseqDotHier.pdf", sep = "-"),width = 15, height = 8)
Clustered_DotPlot(rds, features = abseq, colors_use_exp = rev(colors(20)), 
                  cluster_ident = F, cluster_feature = T,
                  exp_color_min = 0, show_ident_legend = T)
dev.off()
DefaultAssay(rds) <- "SCT"

## Ridge ####
group <-"geno"
Idents(rds) <- group
DefaultAssay(rds) <- "ADT"
abseq <- rownames(rds@assays$ADT)
plots <- RidgePlot(rds, features = abseq, stack = T) + NoLegend()
# PRINT
#pdf(paste(rds@project.name, "AbseqRidge.pdf", sep = "-"), width = length(abseq) / 5 * 4.8, height = length(abseq) / 5 * 6)
pdf(paste0(rds@project.name, "-AbseqRidge-", group, ".pdf"), width = 12, height = 3)
print(plots + labs(title = paste("Abseq Ridge Plot:", group)) + 
        theme(axis.text.y = element_text(hjust = 0))
)
dev.off()

DefaultAssay(rds) <- "SCT"

## Featplot #### 
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

## Violin #############
# PRINT
DefaultAssay(rds) <-"ADT"
pdf(paste0(rds@project.name, "AbVln.pdf"), width = 28, height = 4)
print(VlnPlot(rds, features = abseq[0:8], pt.size = 0, ncol = 8))
print(VlnPlot(rds, features = abseq[9:16], pt.size = 0,  ncol = 8))
print(VlnPlot(rds, features = abseq[17:24], pt.size = 0, ncol = 8))
print(VlnPlot(rds, features = abseq[25:31], pt.size = 0, ncol = 7))
dev.off()
DefaultAssay(rds) <-"SCT"

## Regular Heatmap ####
library(pheatmap)
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
# PRINT
pdf(paste(rds@project.name, "AbsHeat.tiss.geno.pdf", sep = "-"), width = 4, height = 8)
avg_exp <- AggregateExpression(rds, assays = "ADT", slot = "data")
avg_matrix <- avg_exp$ADT[abseq, ]
scaled <- t(scale(t(avg_matrix)))
print(pheatmap(scaled, cluster_rows = T, cluster_cols = T, main = "Abseq Avg Expression"))
dev.off()
DefaultAssay(rds) <- "SCT"

# GENES ######################### 
### CHECK ####
grep("GBP", rownames(rds), value = T)
#write.csv(grep("Il", rds_genes, value = T), file = "cytokinesInBraindata.csv")

## DotPlot ####
listName <- "Type1_2_dress"
group <- "anno.cond"
w = 10
h = 6
n = 2
Idents(rds) <- group
table(Idents(rds))
DefaultAssay(rds) <- "SCT"
genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)

# PRINT
pdf(paste0(rds@project.name , "-Dot-", listName, "-", group, ".pdf"), width = w, height = h)
for (col in names){
  print(col)
  print(DotPlot(rds, features = toupper(genes[[col]]), cols = "RdYlBu", col.min = 0, dot.scale = 5) +
          #coord_flip() + 
          geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
          ggtitle(paste0(col, " | width: ", w, " height: ", h)) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) + 
        # geom_hline(yintercept = seq(n-1.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black", linetype = "dashed") +
        geom_hline(yintercept = seq(n+0.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black")
  )
}
dev.off()

#### scCustom ####
colors <- colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))
listName <- "Microglia-Dotlist-anno.geno_3"
group <- "anno.geno"
Idents(rds) <- group
DefaultAssay(rds) <- "SCT"
genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)

# PRINT
pdf(paste(rds@project.name, listName, group, "DotSCC.pdf", sep = "-"), height =6, width = 15)
for (col in names){
  print(col)
  plot.new()
  text(x = 0.5, y = 0.5, col, cex= 2)
  Clustered_DotPlot(rds, features = unlist(genes[col]), flip = T, 
                    colors_use_exp = rev(colors(20)), x_lab_rotate = T, 
                    group.by = group, exp_color_min = 0, 
                    plot_km_elbow = F)
}
dev.off()

## Ridge ####
listName <- "Efferocytosis"
Idents(rds) <- rds$seurat_clusters
DefaultAssay(rds) <- "SCT"
genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)

# PRINT
pdf(paste(rds@project.name, "Ridge.pdf", sep = "-"),width = 25, height = 6)
for (col in names) {
  print(col)
  plots <- RidgePlot(rds, features = test, stack = T)
  print(plots + plot_annotation((title =col)) + NoLegend())
  
}
dev.off()

## FeatPlot #### 
listName <- "NeutProgenitor"

DefaultAssay(rds) <- "SCT"
genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)

# PRINT
#lapply(genes[col], toupper)
pdf(paste(rds@project.name, listName, "Feature-RNA.pdf", sep = "-"), width = 12, height = 12)
for (col in names){
  print(col)
  plots <- FeaturePlot(rds, reduction = rds@misc$subclustUmap, features = genes[[col]], ncol = 3, raster = F) & 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(), 
          axis.ticks = element_blank())
  plots <- wrap_plots(plots) + plot_annotation(title =col)
  print(plots)
  rm(plots)
}
dev.off()

#### scCustom AddModule #####
# Get genes with prefix
prefix <- "HBB"
groupName <- "RBC"
rds_genes <- rownames(rds)
genes <- grep(prefix, rds_genes, value = T)
genes
rds <- AddModuleScore(rds, features = genes, name = "group")

# When plotting, features must be the group score name with "1" at the end
pdf(paste0(rds@project.name, "-Feature-", prefix, ".pdf"))
FeaturePlot(rds, features = "group1", raster = F) + ggtitle(prefix)
dev.off()

# scCustom density featureplot
library(scCustomize)
library(viridis)
pal <- viridis(n = 10, option = "D")

listName <- "NeutProgenitor"

genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)
names

Idents(rds) <- rds$tissue.geno
idents <- levels(unique(Idents(rds)))
pdf(paste0(rds@project.name, "-Feature-", listName, ".pdf"))
for (ident in idents) {
  sample <- subset(rds, idents = ident)
  for (col in names){
    print(col)
    sample <- AddModuleScore(sample, features = genes[col], name = "group")
    print(FeaturePlot_scCustom(sample, features = "group1", reduction = "rna.umap") + ggtitle(ident))
  }
}
dev.off()

## Coexpression #######
pdf(paste0(rds@project.name, "-Coexpress-Eos.pdf"), width = 15, height = 4)
FeaturePlot(rds, features = c("Ear2", "adt_SiglecF"), reduction = rds@misc$umap, blend = T)
dev.off()

## Regular Heatmap ####
library(pheatmap)

listName <- "Neut-trx"
group <- "anno.geno"
n = 4
Idents(rds) <- rds[[group]][,1]
DefaultAssay(rds) <- "SCT"
genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)

# PRINT
pdf(paste(rds@project.name, listName, group, "Heat.pdf", sep = "-"), width = 7,height = 15)
for (col in names){
  print(col)
  avg_exp <- AggregateExpression(rds, group.by = group, assays = "SCT", slot = "data")
  avg_matrix <- avg_exp$SCT[genes[[col]], ]
  scaled <- t(scale(t(avg_matrix)))
  print(pheatmap(scaled, cluster_rows = F, cluster_cols = F, main = "Average Expression per Cluster"))
}
dev.off()

## Violin ####
DefaultAssay(rds) <- "SCT"
listName <- "Yoon2025-Inflammasomes"

genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>%
  lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)
pdf(paste(rds@project.name, listName, "RNAVln.pdf", sep = "-"), width = 4, height = 6)
for (col in names){
  print(col)
  plot <- VlnPlot(rds, features = genes[[col]], group.by = group, pt.size = 0, stack = T, flip = T)
  print(plot + xlab("") + ylab("") + NoLegend() + ggtitle(col))
}
dev.off()

# DEG ###################################
## Cluster V. All ####
# rds <- PrepSCTFindMarkers(rds)
# print("PrepSCTFindMarkers done")

group <-"clust.geno"
Idents(rds) <- group
idents <- levels(Idents(rds))
idents
degs <- FindAllMarkers(rds, only.pos = F, logfc.threshold = 0.5, assay = "SCT", recorrect_umi = F)
dim(degs)
write.csv(degs, paste0(rds@project.name, "-", group, "-DEGs.csv"))

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


## Filter DEGs ####
# by LogFC, gate out upregs if pct1 < 0.5
### IF BLANK SHEETS MANUALLY DELETE!!!!!!
input_file <- paste0(rds@project.name, "-", group, "-DEGs")

# Also dotplot filtered genes
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

## Dotplot filtered DEGs ####
colors <- colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))
pdf(paste(input_file, "DotTop100.pdf", sep = "-"), width = 20, height = 6)
for (sheet in sheet_names) {
  # Read the sheet
  data <- read_excel(paste0(input_file, "-filtered.xlsx"), sheet = sheet)
  # Print cluster name
  print(sheet)
  plot.new()
  # x and y coordinates are from 0 (left/bottom) to 1 (right/top)
  text(x = 0.5, y = 0.5, paste(sheet, "100 Upreg"), cex = 5, font = 2) # cex changes text size, font changes style
  
  Clustered_DotPlot(rds, features = pull(data, gene), flip = T, 
                    colors_use_exp = rev(colors(20)), x_lab_rotate = T, 
                    group.by = group, exp_color_min = 0, 
                    plot_km_elbow = F)
}
dev.off()

## DEG.xlsx -> Volcano ####
#input_file_list <- c("T Cell-GCSF-IL6-DEGs-wilcox", "Monocyte-GCSF-IL6-DEGs-wilcox")

# for (input_file in input_file_list){
# input_file = "Anno-dressWB-basicAnno-DEGs"
# print(input_file)

sheet_names <- excel_sheets(paste0(input_file, ".xlsx"))
pdf(paste0(input_file, "-Volcano.pdf"))
for (sheet in sheet_names){
  print(sheet)
  degs <- read.xlsx(paste0(input_file, ".xlsx"), sheet = sheet)
  print(EnhancedVolcano(degs, lab = degs$gene, x = "avg_log2FC", y = "p_val_adj",
                        pCutoff = 0.05, FCcutoff = 2.0,
                        title = sheet, subtitle = ""))
}
dev.off()
# }

## Dotplot top 10 DEGs ####
# FILTER OUT SMALL CLUSTERS
# n <- 50
rds$clust.cond <- paste(rds$seurat_clusters, rds$condition)
Idents(rds) <- rds$clust.cond
# clustSize <- table(rds$clust.cond)
# keep  <- names(clustSize[clustSize >= n])
# sub <-subset(rds, idents = keep)
### by cluster ####
top <- 10
group <- "seurat_clusters"
input_file <- paste0(rds@project.name, "-", group,"-DEGs")
sheet_names <- excel_sheets(paste0(input_file, "-filtered.xlsx"))
genes <- c()
for (sheet in sheet_names) {
  print(sheet)
  # Read the sheet
  degs <- read.xlsx(paste0(input_file, "-filtered.xlsx"), sheet = sheet)
  # Filter out noncoding
  degs <- degs[!grepl("^ENS", degs$gene),]
  degs <- degs[!grepl("^LINC", degs$gene),]
  degs <- degs[!grepl("Rik$",degs$gene),]
  genes <- append(genes, degs$gene[1:top])
}
genes <- unique(genes[!is.na(genes)])

pdf(paste0(rds@project.name,"-", group, "-DotTop", top,".pdf"), width =17, height =5)
print(DotPlot(rds, features = unique(genes), cols = "RdYlBu", col.min = 0, dot.scale = 5, group.by = group) +
        ggtitle(paste(rds@project.name, "Top", top, "per", group)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) +
        geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5))
dev.off()


# Pathway ####
## DEenrich ####
# https://satijalab.org/seurat/reference/deenrichrplot
# install.packages("enrichR)
DEenrichRPlot(rds, return.gene.list = T,
              # ident.1 = "Brain Mature-Active Neutrophil",
              # ident.2 = "Brain Active Neutrophil",
              enrich.database = "GO_Biological Process",
              max.genes = 50)

## clusterProfiler ####
# Pathway Enrichment Workflow for Seurat FindAllMarkers Output
# - Loads marker table
# - Builds per-cluster gene lists
# - Runs GO enrichment with clusterProfiler
# - Generates dotplots, barplots, and multi-cluster comparisons
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)   # org.Mm.eg.db for mouse, org.Hs.eg.db human
library(ReactomePA)

## 1. Load FindAllMarkers output
# Expecting columns: gene, cluster, avg_log2FC, p_val_adj, etc.
markers <- read.csv(paste0(rds@project.name,"DEGs.csv"))

## 2. Filter DEGs and split into per-cluster gene lists ------------------------
# Adjust thresholds as needed
gene_lists <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  split(.$cluster) %>%
  lapply(function(df) df$gene)

## 3. Run GO enrichment for each cluster ---------------------------------------
# ont = "BP" (Biological Process), can also use "MF" or "CC"
ego_list <- lapply(gene_lists, function(genes) {
  enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
})

## 4. Plot enrichment for a single cluster -------------------------------------
# Example: cluster "0" (change index as needed)
cluster_id <- names(ego_list)[1]

# Dotplot (most common)
dotplot(ego_list[[cluster_id]], showCategory = 20) +
  ggtitle(paste("GO BP Enrichment – Cluster", cluster_id))

# Barplot
barplot(ego_list[[cluster_id]], showCategory = 20) +
  ggtitle(paste("Top GO Terms – Cluster", cluster_id))

# Enrichment map (clusters similar pathways)
emapplot(pairwise_termsim(ego_list[[cluster_id]]))

# Gene–pathway network
cnetplot(ego_list[[cluster_id]], showCategory = 10)

## 5. Multi-cluster comparison --------------------------------------------------
# Produces a single dotplot comparing pathways across clusters
cc <- compareCluster(
  geneCluster = gene_lists,
  fun         = "enrichGO",
  OrgDb       = org.Hs.eg.db,
  ont         = "BP"
)

dotplot(cc, showCategory = 10) +
  ggtitle("GO BP Enrichment Across Clusters")

## 6. Optional: Reactome or KEGG -----------------------------------------------
# Reactome
react_list <- lapply(gene_lists, function(genes) {
  enrichPathway(gene = genes, organism = "human", pvalueCutoff = 0.05)
})

# KEGG
kegg_list <- lapply(gene_lists, function(genes) {
  enrichKEGG(gene = genes, organism = "hsa", pvalueCutoff = 0.05)
})

# Example Reactome dotplot
dotplot(react_list[[cluster_id]], showCategory = 20) +
  ggtitle(paste("Reactome Pathway Enrichment – Cluster", cluster_id))


# MISC ####
## transfer annotations ####
# Get cell IDs
oldID <- colnames(old)
newID <- colnames(rds)
shared <- intersect(oldID, newID)
length(shared)

# Get annotations
oldAnno <- old@meta.data[shared, "annotations"]
rds@meta.data[shared, "oldAnno"] <- oldAnno
table(rds$oldAnno)

# Print
Idents(rds) <- rds$oldAnno
pdf("Nomid-WNN0.1-oldAnno.pdf",width = 15, height = 8)
DimPlot(rds, reduction = "wnn.umap", label = T, raster = F)
dev.off()


# Sample Group Naming ####################################################
# alphabetize the cell types
Idents(rds) <- rds$anno.orig
Idents(rds) <- factor(Idents(object = rds), levels = sort(levels(rds)))

## NOMID UNST ####
rds$geno <- case_when(
  rds$orig.ident == "C53-nomid-exact-poly10-51x71" ~ "WT",
  rds$orig.ident == "C54-nomid-exact-poly10-51x71" ~ "WT",
  rds$orig.ident == "C55-nomid-exact-poly10-51x71" ~ "NM",
  rds$orig.ident == "C56-nomid-exact-poly10-51x71" ~ "NM",
)

rds$tissue <- case_when(
  rds$orig.ident == "C53-nomid-exact-poly10-51x71" ~ "Brain",
  rds$orig.ident == "C54-nomid-exact-poly10-51x71" ~ "Spleen",
  rds$orig.ident == "C55-nomid-exact-poly10-51x71" ~ "Brain",
  rds$orig.ident == "C56-nomid-exact-poly10-51x71" ~ "Spleen",
)
rds$tiss.geno <- paste(rds$tissue, rds$geno)
table(rds$geno)
table(rds$tissue)
table(rds$tiss.geno)

## NOMID #####
rds$orig.ident <- recode(rds$orig.ident, "NOMID-SPLEEN" = "spleen-NM", "WT-SPLEEN" = "spleen-WT")
rds$anno.orig = paste0(rds$annotations, "_", rds$orig.ident)

rds$geno <- case_when(
  rds$orig.ident == "C53-nomid-exact-poly10-51x71" ~ "WT",
  rds$orig.ident == "C54-nomid-exact-poly10-51x71" ~ "WT",
  rds$orig.ident == "C55-nomid-exact-poly10-51x71" ~ "NM",
  rds$orig.ident == "C56-nomid-exact-poly10-51x71" ~ "NM",
  rds$orig.ident == "C77-nomid-exact-poly10-51x71" ~ "WT",
  rds$orig.ident == "C78-exact-poly10-51x71" ~ "WT",
  rds$orig.ident == "C79-exact-poly10-51x71" ~ "WT",
  rds$orig.ident == "C80-exact-poly10-51x71" ~ "WT",
  rds$orig.ident == "C81-exact-poly10-51x71" ~ "NM",
  rds$orig.ident == "C82-exact-poly10-51x71" ~ "NM",
  rds$orig.ident == "C83-nomid-exact-poly10-51x71" ~ "NM",
  rds$orig.ident == "C84-exact-poly10-51x71" ~ "NM",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)
table(rds$geno)

rds$mouse <- case_when(
  rds$orig.ident == "C53-nomid-exact-poly10-51x71" ~ "WT",
  rds$orig.ident == "C54-nomid-exact-poly10-51x71" ~ "WT",
  rds$orig.ident == "C55-nomid-exact-poly10-51x71" ~ "NM",
  rds$orig.ident == "C56-nomid-exact-poly10-51x71" ~ "NM",
  rds$orig.ident == "C77-nomid-exact-poly10-51x71" ~ "WT1",
  rds$orig.ident == "C78-exact-poly10-51x71" ~ "WT2",
  rds$orig.ident == "C79-exact-poly10-51x71" ~ "WT3",
  rds$orig.ident == "C80-exact-poly10-51x71" ~ "WT4",
  rds$orig.ident == "C81-exact-poly10-51x71" ~ "NM1",
  rds$orig.ident == "C82-exact-poly10-51x71" ~ "NM2",
  rds$orig.ident == "C83-nomid-exact-poly10-51x71" ~ "NM3",
  rds$orig.ident == "C84-exact-poly10-51x71" ~ "NM4",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)
table(rds$mouse)

rds$tissue <- case_when(
  rds$orig.ident == "C53-nomid-exact-poly10-51x71" ~ "Brain",
  rds$orig.ident == "C54-nomid-exact-poly10-51x71" ~ "Spleen",
  rds$orig.ident == "C55-nomid-exact-poly10-51x71" ~ "Brain",
  rds$orig.ident == "C56-nomid-exact-poly10-51x71" ~ "Spleen",
  rds$orig.ident == "C77-nomid-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag04_mm" ~ "Brain",
  rds$orig.ident == "C77-nomid-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag05_mm" ~ "Spleen",
  rds$orig.ident == "C78-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag04_mm" ~ "Brain",
  rds$orig.ident == "C78-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag05_mm" ~ "Spleen",
  rds$orig.ident == "C79-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag04_mm" ~ "Brain",
  rds$orig.ident == "C79-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag05_mm" ~ "Spleen",
  rds$orig.ident == "C80-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag06_mm" ~ "Brain",
  rds$orig.ident == "C80-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag07_mm" ~ "Spleen",
  rds$orig.ident == "C81-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag06_mm" ~ "Brain",
  rds$orig.ident == "C81-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag07_mm" ~ "Spleen",
  rds$orig.ident == "C82-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag06_mm" ~ "Brain",
  rds$orig.ident == "C82-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag07_mm" ~ "Spleen",
  rds$orig.ident == "C83-nomid-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag08_mm" ~ "Brain",
  rds$orig.ident == "C83-nomid-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag09_mm" ~ "Spleen",
  rds$orig.ident == "C84-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag08_mm" ~ "Brain",
  rds$orig.ident == "C84-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag09_mm" ~ "Spleen",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)
table(rds$tissue)

rds$tissue.geno = paste0(rds$tissue, " ", rds$geno)
rds$tissue.mouse = paste0(rds$tissue, " ", rds$mouse)
table(rds$tissue.geno)
table(rds$tissue.mouse)


## DRESS WB ####
rds$patient <- case_when(
  rds$orig.ident == "C36-exact-poly10-51x71" ~ "sick-GGW",
  rds$orig.ident == "C47-1-exact-poly10-51x71" ~ "reco-GGW",
  rds$orig.ident == "C47-2-exact-poly10-51x71" ~ "reco-GGW",
  rds$orig.ident == "C85-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag02_hs" ~ "preDupiAsth221",
  rds$orig.ident == "C85-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag03_hs" ~ "sick-BS",
  rds$orig.ident == "C86-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag04_hs" ~ "sick-PDJ-BF",
  rds$orig.ident == "C86-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag05_hs" ~ "sick-PDJ",
  rds$orig.ident == "C92-exact-poly10-51x71" ~ "reco-PDJ",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)

rds$condition <- case_when(
  rds$orig.ident == "C36-exact-poly10-51x71" ~ "Dress",
  rds$orig.ident == "C47-1-exact-poly10-51x71" ~ "recovered",
  rds$orig.ident == "C47-2-exact-poly10-51x71" ~ "recovered",
  rds$orig.ident == "C85-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag03_hs" ~ "Dress",
  rds$orig.ident == "C86-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag04_hs" ~ "Dress-BF",
  rds$orig.ident == "C86-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag05_hs" ~ "Dress",
  rds$orig.ident == "C92-exact-poly10-51x71" ~ "recovered",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)
table(rds$patient)
table(rds$condition)

## DRESS PBMC ####
rds$patient <- case_when(
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_hs" ~ "DRS01",
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_hs" ~ "DRS01",
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_hs" ~ "DRS01",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag08_hs" ~ "DRS02",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag09_hs" ~ "DRS02",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_hs" ~ "DRS03",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_hs" ~ "DRS03",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_hs" ~ "DRS03",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)
table(rds$patient)

rds$condition <- case_when(
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_hs" ~ "active",
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_hs" ~ "treated",
  rds$orig.ident == "C104-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_hs" ~ "recovered",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag08_hs" ~ "active",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag09_hs" ~ "recovered",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_hs" ~ "active",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_hs" ~ "treated",
  rds$orig.ident == "C105-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_hs" ~ "recovered",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)
table(rds$condition)

## COCCI ####
rds$strain.group = paste0(rds$strain, "_", rds$group )
rds$clust.group.strain = paste0(rds$seurat_clusters, "_", rds$group, "_", rds$strain)
Idents(rds) <- rds$clust.group.strain

# Enrichr formatting ####
library(stringr)

enrichrFormat <- function(text){
  sep <- unlist(strsplit(text, ";"))
  for (gene in sep){
    write.csv(str_to_title(sep), "genes.csv")
  }
}

enrichrFormat("C1QB;CX3CR1;P2RY12;C1QA;CSF1R;P2RY13;CD83;GPR34;AIF1;HPGDS;SALL1;OLFML3;CCL4;MRC1;C3AR1;TMEM119;SLC29A3;C1QC")
