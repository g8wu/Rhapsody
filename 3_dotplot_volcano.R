#install.packages('devtools')
devtools::install_github('immunogenomics/presto')
BiocManager::install('multtest')
install.packages('metap')
set.seed(99)
if(!require("EnhancedVolcano", quietly = T)) BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
if(!require("ComplexHeatmap", quietly = T)) BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(Seurat)
library(patchwork)
library(ggplot2)
if(!require("scCustomize", quietly = T)) install.packages("scCustomize")
library(scCustomize)
library(RColorBrewer)
if(!require("ggpubr", quietly = T)) install.packages("ggpubr")
library(ggpubr)
library(dplyr)
display.brewer.all(colorblindFriendly = TRUE)

res <- "0.3"
project <- "SPLEEN"
fileType <- "nomid-exact-p10-51x71"

#############  QC post Cluster ############# 
features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_ADT")

pdf(paste(project, fileType, res, "RNAClustVlnPlot.pdf", sep = "-"), width = 20, height = 6)
print(VlnPlot(rds, features = features, pt.size = 0, ncol = length(features), group.by = "orig.ident"))
print(VlnPlot(rds, features = features, pt.size = 0, ncol = length(features), group.by = "seurat_clusters"))
dev.off()

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

pdf(paste(project, fileType, res, "UMIBarplot.pdf", sep = "-"))
print(umiBar)
dev.off()

#############  Abseq Dotplot ############# 
abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
colors <- colorRampPalette(brewer.pal(n = 9, name = "Oranges"))
pdf(paste(project, fileType, res, "AbseqDotHierOrange.pdf", sep = "-"),width = 15, height = 10)
Clustered_DotPlot(rds, features = abseq, colors_use_exp = colors(20), 
                  exp_color_min = 0, show_ident_legend = F)
dev.off()
DefaultAssay(rds) <- "SCT"

############# Abseq Ridge Plots
# PRINT
plots <- RidgePlot(rds, features = abseq, ncol = 5)
# PRINT
pdf(paste(project, fileType, res, "AbseqRidge.pdf", sep = "-"), width = length(abseq) / 5 * 4.8, height = length(abseq) / 5 * 5)
print(plots)
dev.off()

############# Abseq Featureplots ############# 
plots <- FeaturePlot(rds, reduction = "rna.umap", features = abseq, ncol = 5) & 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank())
# PRINT
pdf(paste(project, fileType, res, "AbseqFeaturePlots.pdf", sep = "-"), width = length(abseq) / 5 * 4.8, height = length(abseq) / 5 * 5)
print(plots)
dev.off()

############# GENES Dot Plot ############# 
rds$orig.ident <- recode(rds$orig.ident, "NOMID-SPLEEN" = "spleen-NM", "WT-SPLEEN" = "spleen-WT")
rds$orig.clust = paste0(rds$seurat_clusters, "_", rds$orig.ident )

# Checking for genes with different names
rds_genes <- rownames(GetAssayData(rds, assay = "RNA"))
grep("S100a8", rds_genes, value = T)
write.csv(grep("Il", rds_genes, value = T), file = "cytokinesInBraindata.csv")

# START!!
DefaultAssay(rds) <- "SCT"
genes <- read.csv("CellType.csv", header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)
pdf(paste(project, fileType, res, "CellTypeList-SCTDot.pdf", sep = "-"))
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
pdf(paste(project, fileType, res, "InterestList-SCTDot.pdf", sep = "-"), height = 11, width = 13)
for (col in names){
  print(col)
  print(DotPlot(rds, features = genes[[col]], col.min = 0, group.by = "orig.clust") + 
          labs(title = paste("SCT Dotplot:", col)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}
dev.off()

####### DEGs & Volcanoes for each cluster #######
rds <- PrepSCTFindMarkers(rds)
saveRDS(rds, file = paste(project, fileType, res, "DEGs.RDS", sep = "-"))
degs <- FindAllMarkers(rds, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(degs, file=paste(project, fileType, res, "DEGs.csv", sep = "-"))

pdf(paste(project, fileType, res, "Volcano.pdf", sep = "-"))
#for (i in 1:length(unique(rds$seurat_clusters))){
for (i in 1:2){
  clustDegs <- subset(degs, cluster == i)
  print(EnhancedVolcano(clustDegs, rownames(clustDegs), x = "avg_log2FC", y = "p_val_adj",
                        title = paste(project, "cluster", i, sep = " "), subtitle = ""))
}
dev.off()

# degs <- degs[order(degs$avg_log2FC, decreasing = TRUE), ]
# topDegs <- rbind(head(degs, 100), tail(markers, 100))
# write.csv(topDegs, file=paste(project, fileType, res,"DEGs.csv", sep = "-"))
# 
# #how cluster 9 changes wrt cluster 12
# rds <- FindMarkers(rds, ident.1 = "9", ident.2 = "12", min.pct=0.25, logfc.threshold=0.1) 
# 
# clusters <- unique(rds$seurat_clusters)
# clusters <- c(1)
# for (i in clusters){
#   print(i)
#   DefaultAssay(rds) <- "SCT"
#   #rds <- PrepSCTFindMarkers(rds)
#   markers <- FindConservedMarkers(rds, assay = "SCT", ident.1 = i, grouping.var = "orig.ident", verbose = FALSE)
#   markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]
#   markers <- markers[markers$p_val < 0.05, ]
#   markers <- rbind(head(markers, 100), tail(markers, 100))
#   write.csv(markers, file=paste(project, fileType, res, i, "DEGs.csv", sep = "-"))
#   pdf(paste(project, fileType, res, "Volcano", i, ".pdf", sep = "-"), height = 10, width = 10)
#   print(EnhancedVolcano(markers, rownames(markers), x = "avg_log2FC", y = "p_val_adj",
#                         title = paste(project, i, sep = " "), subtitle = ""))
#   dev.off()
# }
# 
# pdf(paste("Brain", i,"Volcanoplots.pdf", sep = "-"), height = 10, width = 10)
# print(EnhancedVolcano(degs, rownames(degs), x = "avg_log2FC", y = "p_val_adj",
#                       title = paste(project), subtitle = ""))
# dev.off()
# print(paste("brain", i, "done"))
# rm(temp)
# gc()


####### GENES Violin Plot #######
# # START!!
# DefaultAssay(rds) <- "RNA"
# genes <- read.csv("CellTypeList.csv", header = T, na.strings = "") %>%
#   lapply(function(column) {column[!is.na(column) & column != ""]})
# names <- names(genes)
# for (col in names){
#   pdf(paste(project, fileType, res, "RNAVln.pdf", sep = "-"), width = 8.5, height = 11)
#   print(col)
#   plot <- VlnPlot(rds, features = genes[[col]], group.by = "seurat_clusters", pt.size = 0, stack = T, flip = T)
#   print(plot + xlab("") + ylab("") + NoLegend() + ggtitle(col))
#   dev.off()
# }

####### GENES Heat Map #######
# DefaultAssay(rds) <- "SCT"
# genes <- read.csv("InterestList.csv", header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
# names <- names(genes)
# 
# pdf(paste(project, fileType, res, "InterestGenes-SCTHeat.pdf", sep = "-"))
# for (col in names){
#   print(col)
#   print(DoHeatmap(rds, features = genes[[col]], size = 3))
# }
# dev.off()

############# Abseq Violin Plots
# PRINT
# pdf(paste0(project, fileType, res, "AbseqViolin.pdf", sep = "-"), width = 28, height = 4)
# print(VlnPlot(rds, features = abseq[0:8], pt.size = 0, ncol = 8))
# print(VlnPlot(rds, features = abseq[9:16], pt.size = 0,  ncol = 8))
# print(VlnPlot(rds, features = abseq[17:24], pt.size = 0, ncol = 8))
# print(VlnPlot(rds, features = abseq[25:31], pt.size = 0, ncol = 7))
# dev.off()

####### mito/cocci Featureplots #######
# pdf(paste0(project, res, "mitoReads.pdf"))
# print(FeaturePlot(rds, features = "percent.mt", reduction = "rna.umap"))
# dev.off()
# 
# pdf(paste(project, fileType, res, "cocciReads.pdf", sep = "-"))
# print(FeaturePlot(rds, features = "percent.cocci", reduction = "rna.umap"))
# dev.off()

