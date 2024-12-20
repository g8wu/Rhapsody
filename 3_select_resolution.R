res <- "0.5"
project <- "BALF-Pilot"
fileType <- "mouseOnly-poly10-51x71"
wkdir <- getdir()
rds <- readRDS(paste0(project, "-", fileType, "-WNN", res, ".rds"))


write.csv(t(table(Idents(rds))), file=paste0(wkdir, "/", project, "-", res, "_clusterTable.csv"))
pdf(paste0(project, "-", fileType, res, "_bySampleUMAP.pdf"))
print(DimPlot(rds, reduction = 'umap', label.size = 2.5) +
        NoLegend() + labs(title = paste("UMAP Louvian res:", res, ":", project)))
idents <- unique(rds$orig.ident)
for (i in seq_along(idents)){
  sample <- subset(rds, subset = orig.ident == idents[i])
  print(DimPlot(rds, reduction = "wnn.umap", cells.highlight = Cells(sample)) + NoLegend() + 
          labs(title = paste("WNN UMAP Louvian res:", res, ":", project, " ", idents[i])))
}
dev.off()

#############  Abseq Dotplot
DefaultAssay(rds) <- "SCT"
abseq <- rownames(rds[["ADT"]])
# PRINT
pdf(paste0(project,"-",fileType, "-", res, "_AbseqDotplot.pdf"), width = 28, height = 11)
print(DotPlot(rds, features = abseq, cols = "RdBu", col.min = -1, dot.scale = 8, assay = "ADT") + 
        RotatedAxis() + labs(title = paste("AbSeq Dotplot:", project)))
dev.off()
DefaultAssay(rds) <- "RNA"

############# Abseq Featureplots
# PRINT
pdf(paste0(project, "-", fileType, res, "_AbseqFeaturePlots.pdf"), width = 28, height = 4)
print(FeaturePlot(rds, features = abseq[0:8], reduction = 'rna.umap', ncol = 8) & NoLegend())
print(FeaturePlot(rds, features = abseq[9:16], reduction = 'rna.umap',  ncol = 8) & NoLegend())
print(FeaturePlot(rds, features = abseq[17:24],  reduction = 'rna.umap', ncol = 8) & NoLegend())
print(FeaturePlot(rds, features = abseq[25:31],  reduction = 'rna.umap', ncol = 7) & NoLegend())
dev.off()
write(paste0(Sys.time(), " -> READY FOR ANNOTATION!\n"), file = LOG, append = TRUE)

############# Abseq Violin Plots
# PRINT
pdf(paste0(project, "-", fileType, res, "_AbseqViolin.pdf"), width = 28, height = 4)
print(VlnPlot(rds, features = abseq[0:8], pt.size = 0, ncol = 8))
print(VlnPlot(rds, features = abseq[9:16], pt.size = 0,  ncol = 8))
print(VlnPlot(rds, features = abseq[17:24], pt.size = 0, ncol = 8))
print(VlnPlot(rds, features = abseq[25:31], pt.size = 0, ncol = 7))
dev.off()

############# Abseq Ridge Plots
# PRINT
pdf(paste0(project, "-", fileType, res, "_AbseqRidge.pdf"), width = 28, height = 4)
print(RidgePlot(rds, features = abseq[0:8], ncol = 8))
print(RidgePlot(rds, features = abseq[9:16],  ncol = 8))
print(RidgePlot(rds, features = abseq[17:24], ncol = 8))
print(RidgePlot(rds, features = abseq[25:31], ncol = 7))
dev.off()

####### DEGs for each cluster #######
n <- length(unique(Idents(rds)))-1
degs <- list()
for (i in 0:n){
  markers <- FindMarkers(rds, ident.1 = i, min.pct = 0.15, verbose = FALSE)
  markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]
  degs <- rbind(head(markers, 100), tail(markers, 100))
  write.csv(degs, file=paste0("Cluster_", as.character(i), "_DEGs.csv"))
  print(i)
}  

####### GENES Violin Plot #######
# Checking for genes with different names
rds_genes <- rownames(GetAssayData(rds, assay = "RNA"))
grep("Fgcr", rds_genes, value = T)

# START!!
genes <- read.csv("CellTypeList.csv", header = T, na.strings = "") %>%
  lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)
pdf(paste0(project, "CellTypeViolin.pdf"))
for (col in names){
  print(col)
  plots <- lapply(genes[[col]], function(gene) VlnPlot(rds, features = gene, group.by = "seurat_clusters") + 
                    NoLegend() + xlab("") + ylab(""))
  combined_plots <- wrap_plots(plots, ncol = 2) + plot_annotation(title = col, theme = theme(plot.title = element_text(face = "bold")))
  print(combined_plots)
}
dev.off()

genes <- read.csv("InterestList.csv", header = T, na.strings = "") %>%
  lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)
pdf(paste0(project, "_InterestGenesViolin.pdf"))
for (col in names){
  print(col)
  plots <- lapply(genes[[col]], function(gene) VlnPlot(rds, features = gene, group.by = "seurat_clusters") + 
                    NoLegend() + xlab("") + ylab(""))
  combined_plots <- wrap_plots(plots, ncol = 2) + plot_annotation(title = col, theme = theme(plot.title = element_text(face = "bold")))
  print(combined_plots)
}
dev.off()