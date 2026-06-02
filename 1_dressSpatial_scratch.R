set.seed(99)
# VlnPlot ####
group <- "annotations"
Idents(rds) <- group

pdf(paste(project, "VlnQC.pdf", sep = "-"), width = 20, height = 7)
print(VlnPlot(rds, features = features, pt.size = 0, ncol = 4, group.by = "patient") + NoLegend())
print(VlnPlot(rds, features = features, pt.size = 0, ncol = 4, group.by = "condition") + NoLegend())
dev.off()

# Norm & Cluster ####
alg <- 4
resRange <- seq(0.2, 0.5, by = 0.1)
rnaPCs <- 50
adtPCs <- 20
algKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")

## RNA Norm PCA ############## 
DefaultAssay(rds) <- 'RNA'
rds <- SCTransform(rds, assay = "RNA",
                   new.assay.name = "SCT",
                   vars.to.regress = "percent.mt")
rds <- RunPCA(rds, assay = "SCT", reduction.name = "pca_rna")

## Harmonize RNA ####
# Batch effect correct (BeC)
rds <- RunHarmony(rds, group.by.vars = "slide", 
                  reduction = "pca_rna", 
                  reduction.save = "harmony_rna")
## UMAP RNA ####
# UMAP Without Harmony
rds <- RunUMAP(rds, reduction = "pca_rna", dims = 1:rnaPCs, reduction.name = "umap")
# UMAP With Harmony
rds <- RunUMAP(rds, reduction = "harmony_rna", dims = 1:rnaPCs, reduction.name = "harmony_umap")

# Print
pdf(paste0(rds@project.name,"-PCA-Harmony.pdf"))
print(DimPlot(rds, reduction = "umap"))
print(DimPlot(rds, reduction = "harmony_umap"))
print(DimPlot(rds, reduction = "pca_rna"))
dev.off()

## RNA Cluster ####
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

## SAVE!! ####
saveRDS(rds, paste0(rds@project.name,".rds"))

# DotPlot ####
listName <- "Type1_2_dress"
rds$anno.cond <- paste(rds$annotations, rds$condition)
Idents(rds) <- "anno.cond"
# alphabetize the cell types
Idents(rds) <- factor(Idents(rds), levels = sort(levels(rds)))
rds$anno.cond <- Idents(rds)

group <- "anno.cond"
n = 2
Idents(rds) <- group
DefaultAssay(rds) <- "SCT"
genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)
# custom <- c(2.5, 4.5, 6.5, 8.5, 10.5, 11.5, 13.5, 15.5, 17.6, 19.5, 21.5, 23.5, 25.5, 27.5)

# PRINT
pdf(paste(rds@project.name, listName, group, "Dot.pdf", sep = "-"), width = 11, height = 8)
for (col in names){
  print(col)
  print(DotPlot(rds, features = toupper(genes[[col]]), cols = "RdYlBu", col.min = 0, dot.scale = 5) +
          #coord_flip() + 
          geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
          ggtitle(col) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) +
          # geom_hline(yintercept = seq(n-1.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black", linetype = "dashed") +
          geom_hline(yintercept = seq(n+0.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black")
  )
}
dev.off()

### Dot celltype ####
listName <- "Mafe_celltype"
group <- "annotations"
Idents(rds) <- group
DefaultAssay(rds) <- "SCT"
genes <- read.csv(paste0(listName, ".csv"), header = T, na.strings = "") %>% lapply(function(column) {column[!is.na(column) & column != ""]})
names <- names(genes)
# custom <- c(2.5, 4.5, 6.5, 8.5, 10.5, 11.5, 13.5, 15.5, 17.6, 19.5, 21.5, 23.5, 25.5, 27.5)

# PRINT
pdf(paste(rds@project.name, listName, group, "Dot.pdf", sep = "-"), width = 11, height = 8)
for (col in names){
  print(col)
  print(DotPlot(rds, features = toupper(genes[[col]]), cols = "RdYlBu", col.min = 0, dot.scale = 5) +
          geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5) +
          ggtitle(col) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0))
        # geom_hline(yintercept = seq(n-1.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black", linetype = "dashed") +
        # geom_hline(yintercept = seq(n+0.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black")
  )
}
dev.off()

# Feature scCustom AddModule #####
rds <- AddModuleScore(rds, features = toupper(genes$IFN.inducible), name = "group")
FeaturePlot(rds, reduction = "SPATIAL", features = "IL17A")

# Mystery clusters ####
> Idents(rds)<-rds$annotations
> clust2 <- WhichCells(rds, idents =2)
> clust4 <- WhichCells(rds, idents =4)
> clust8 <- WhichCells(rds, idents =8)
> clust11 <- WhichCells(rds, idents =11)
> clust13 <- WhichCells(rds, idents =15)
> clust13 <- WhichCells(rds, idents =13)
> clust15 <- WhichCells(rds, idents =15)
> clust16 <- WhichCells(rds, idents =16)
> pdf("mystery-clust-spatial.pdf")
> DimPlot(dressSpatial, reduction = "SPATIAL", cells.highlight = clust2) + ggtitle("Cluster 2")
> dev.off()
null device 
1 
> pdf("mystery-clust-spatial.pdf")
> DimPlot(dressSpatial, reduction = "SPATIAL", cells.highlight = clust2) + ggtitle("Cluster 2")
> ADSpatial <- readRDS("~/dress/spatial/AD_003_09_ConfPositioned_seurat_spatial_merged.rds")
> DimPlot(ADSpatial, reduction = "SPATIAL", cells.highlight = clust4) + ggtitle("Cluster 4")
> DimPlot(ADSpatial, reduction = "SPATIAL", cells.highlight = clust8) + ggtitle("Cluster 8")
> DimPlot(ADSpatial, reduction = "SPATIAL", cells.highlight = clust11) + ggtitle("Cluster 11")
> DimPlot(dressSpatial, reduction = "SPATIAL", cells.highlight = clust13) + ggtitle("Cluster 13")
> DimPlot(ADSpatial, reduction = "SPATIAL", cells.highlight = clust13) + ggtitle("Cluster 13")
> DimPlot(dressSpatial, reduction = "SPATIAL", cells.highlight = clust15) + ggtitle("Cluster 15")
> DimPlot(ADSpatial, reduction = "SPATIAL", cells.highlight = clust15) + ggtitle("Cluster 15")
> DimPlot(dressSpatial, reduction = "SPATIAL", cells.highlight = clust16) + ggtitle("Cluster 16")
> DimPlot(ADSpatial, reduction = "SPATIAL", cells.highlight = clust16) + ggtitle("Cluster 16")
> dev.off()
null device 

#cc <- toupper(genes$CC.Chemokines)

# When plotting, features must be the group score name with "1" at the end
pdf(paste0(rds@project.name, "-Feature-", prefix, ".pdf"))
FeaturePlot(rds,reduction = "SPATIAL", features = "group1", raster = F) #+ ggtitle("Chemokines")
dev.off()

pdf(paste(rds@project.name, "Feats-IFN.pdf", sep = "-"), width = length(genes$Inflammatory.pathways) / 5 * 8, height = length(genes$Inflammatory.pathways) / 5 * 3)
plots <- FeaturePlot(rds, reduction = rds@misc$umap, features = toupper(genes$Inflammatory.pathways), ncol = 5) & 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank())
print(plots+ plot_annotation(title ="Inflammatory Pathways", theme = theme(plot.title = element_text(size =50))))
dev.off()

# Overlay image ####
install.packages("tiff")
library(tiff)

img1 <- readTIFF("/sbgenomics/project-files/1_Analysis/DRESS Research Final-Reid Oldenburg.tif")
img2 <- readTIFF("/sbgenomics/project-files/1_Analysis/DRESS Research Final-Reid Oldenburg (5, x=2863, y=5099, w=43878, h=41463).tif")

## Mirror flip tif ####
if(!require("magick", quietly = T)) install.packages("magick")
library(magick)

im <- image_read("/sbgenomics/project-files/1_Analysis/DRESS Research Final-Reid Oldenburg (5, x=2863, y=5099, w=43878, h=41463).tif")
im_mirror <- image_flop(im)   # horizontal flip
img2 <- im_mirror

## Plot on H&E ####
# H&E in png format
library(png)
library(grid)
library(ggplot2)
library(magick)

df <- cbind(
  as.data.frame(Embeddings(rds, "SPATIAL")),
  annotation = rds$annotations
)
colnames(df)[1:2] <- c("SPATIAL_1", "SPATIAL_2")

img <- readPNG("DRESSslide downsample.png")
g <- rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))
xrange <- range(df$SPATIAL_1)
yrange <- range(df$SPATIAL_2)

### Adjust ####
x_offset <- 200
y_offset <- -250
scale <- 1 # will increase png size by 10%

xmin <- mean(xrange) - diff(xrange)*scale/2 + x_offset
xmax <- mean(xrange) + diff(xrange)*scale/2 + x_offset
ymin <- mean(yrange) - diff(yrange)*scale/2 + y_offset
ymax <- mean(yrange) + diff(yrange)*scale/2 + y_offset

png("HnE1.png", width = 15,height = 15, units = "in", res  = 600)
ggplot(df, aes(SPATIAL_1, SPATIAL_2)) +
  annotation_custom(g, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax) +
  geom_point(aes(color = annotation), size = 0.1, alpha = 1) +
  scale_color_manual(values = scales::hue_pal()(length(unique(df$annotation)))) +
  theme_void() +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 10))) # Big circles in legend
dev.off()

### Rotate ####
img_rot <- image_rotate(image_read("background.png"), 15)  # rotate 15 degrees
g <- rasterGrob(as.raster(img_rot), width = unit(1, "npc"), height = unit(1, "npc"))


# Convert slideseq into seurat ####
counts <- rds@assays$RNA@counts
meta <- rds@meta.data
coords <- rds@images$slice1@coordinates
new_rds <- CreateSeuratObject(counts = counts,
                              meta.data= meta,
                              project = "Dress_FLEX"
)

new_rds@meta.data <- cbind(new_rds@meta.data, coords)

new_rds@misc <- list()
new_rds@tools <- list()
table(new_rds$condition)
saveRDS(new_rds, "DRESS_Flex.rds")

# Combining objects ####
Idents(dressSpatial) <- dressSpatial$condition
rds <- merge(dressSpatial, y = list(AD_Flex_Healthy))
table(Idents(rds))
rds@project.name <- "Dress_healthy_AD"
saveRDS(rds, "dress_healthy.rds")

### Reumap ####
alg <- 4
resRange <- seq(0.2, 0.5, by = 0.1)

# Global variables
rnaPCs <- 50
algKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")

# For laptops/PCs override default variable size limit
n <- 50     # Let variables reach up to n GB
options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB

print("____________________________________________ SCTransform start")
DefaultAssay(rds) <- 'RNA'
rds <- SCTransform(rds, assay = "RNA",
                   new.assay.name = "SCT",
                   vars.to.regress = "percent.mt")
print(paste(Sys.time(), "RNA SCT norm done"))
gc()
### PCA ####
rds <- RunPCA(rds, assay = "SCT", reduction.name = "pca_rna")

### Harmonize RNA ####
# Batch effect correct (BeC)
rds <- RunHarmony(rds, group.by.vars = "slide", 
                  reduction = "pca_rna", 
                  reduction.save = "harmony_rna")
### UMAP RNA ####
# UMAP Without Harmony
rds <- RunUMAP(rds, reduction = "pca_rna", dims = 1:rnaPCs, reduction.name = "umap")
# UMAP With Harmony
rds <- RunUMAP(rds, reduction = "harmony_rna", dims = 1:rnaPCs, reduction.name = "harmony_umap")

# Print
Idents(rds) <- rds$slide
pdf(paste0(rds@project.name,"-PCA-Harmony.pdf"))
print(DimPlot(rds, reduction = "umap"))
print(DimPlot(rds, reduction = "harmony_umap"))
print(DimPlot(rds, reduction = "pca_rna"))
dev.off()

### Multi-res Clust ####
set.seed(99)
rds <- FindNeighbors(rds, dims = 1:rnaPCs, reduction = "pca_rna")

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

# Post Cluster ####
features = c("nFeature_RNA", "nCount_RNA", "percent.mt")

## QC FeatPlot ####
pdf(paste0(rds@project.name, "-QCUMAP.pdf"))
print(FeaturePlot(rds, reduction = "umap", features = features, ncol = 2))
dev.off()

saveRDS(rds, "dress_healthy-BeC.rds")

## Select res ####
rds@misc$umap <- "harmony_umap"
rds$seurat_clusters <- rds$SCT_snn_res.0.2
Idents(rds) <- "seurat_clusters"

# UMAP ####
Idents(rds) <- rds$annotations
pdf(paste0(rds@project.name,"-anno-UMAP.pdf"), width = 7, height =6)
DimPlot(rds, reduction =rds@misc$umap, label = T)
dev.off()

## Dotplot top 10 DEGs ####
rds$clust.cond <- paste(rds$seurat_clusters, rds$condition)
Idents(rds) <- rds$clust.cond

### by cluster ####
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
pdf(paste0(rds@project.name, "-seurat_clusters-DotTop10.pdf"), width = 25, height = 5)
print(DotPlot(rds, features = unique(genes), cols = "RdYlBu", col.min = 0, dot.scale = 5, group.by = "annotations") +
        ggtitle(paste(rds@project.name, "Top 10 per cluster")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y = element_text(hjust = 0)) +
        geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5))
dev.off()


# Subcluster ####
# For laptops/PCs override default variable size limit
n <- 50     # Let variables reach up to n GB
options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB

wkdir <- getwd()
original <- rds
Idents(original) <- original$annotations
idents <- levels(Idents(original))
idents
resRange <- seq(0.05, 0.2, by = 0.05)
alg <- 4
rnaPCs <- 30
algKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")

databaseList <- c("KEGG_2019_Human", 
                  "Jensen_DISEASES_Curated_2025",
                  "Jensen_TISSUES",
                  "MSigDB_Hallmark_2020",
                  "WikiPathways_2024_Human",
                  "Disease_Perturbations_from_GEO_up",
                  "Reactome_Pathways_2024")

# Split cells from cluster back to separate objects by cartridge, rerun integration
for (cluster in idents) {
  set.seed(99)
  # Make new folder for each subcluster
  print(cluster)
  dir.create(paste0(wkdir, "/", cluster))
  setwd(paste0(wkdir,"/", cluster))
  rds <- subset(original, idents = cluster)
  
  ### RNA Normalization #####
  print("____________________________________________ SCTransform start")
  DefaultAssay(rds) <- 'RNA'
  rds <- SCTransform(rds, assay = "RNA",
                     new.assay.name = "SCT",
                     vars.to.regress = "percent.mt")
  print(paste(Sys.time(), "RNA SCT norm done"))
  gc()
  
  ## PCA ####
  rds <- RunPCA(rds, assay = "SCT", reduction.name = "pca_rna")
  
  ## Harmonize RNA ####
  # Batch effect correct (BeC)
  rds <- RunHarmony(rds, group.by.vars = "condition",
                    reduction = "pca_rna",
                    reduction.save = "harmony_rna")
  ### UMAP RNA ####
  # UMAP Without Harmony
  rds <- RunUMAP(rds, reduction = "pca_rna", dims = 1:rnaPCs, reduction.name = "umap")
  # UMAP With Harmony
  rds <- RunUMAP(rds, reduction = "harmony_rna", dims = 1:rnaPCs, reduction.name = "harmony_umap")
  
  # Print
  pdf(paste0(rds@project.name,"-PCA.pdf"))
  Idents(rds) <- "condition"
  print(DimPlot(rds, reduction = "umap"))
  print(DimPlot(rds, reduction = "harmony_umap"))
  print(DimPlot(rds, reduction = "pca_rna"))
  dev.off()
  
  ## Multi-res Clust ####
  ### RNA ONLY ####
  print("_____________________ Starting RNA ONLY Clustering")
  set.seed(99)
  rds <- FindNeighbors(rds, dims = 1:rnaPCs, reduction = "pca_rna")
  
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
  saveRDS(rds, paste0(cluster, ".rds"))
  
  ## QC FeatPlot ####
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
  pdf(paste0(rds@project.name, "-QCUMAP.pdf"))
  print(FeaturePlot(rds, reduction = "harmony_umap", features = features, ncol = 2))
  dev.off()
  
  ## DEG Conditions ####
  rds <- PrepSCTFindMarkers(rds)
  Idents(rds) <- rds$condition
  comparisons <- list(list("DRESS", "Healthy"))
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
  saveWorkbook(wb,file=paste0(rds@project.name,"-condition-DEGs.xlsx"), overwrite = TRUE)
  
  ### Filter DEGs ####
  ### IF BLANK SHEETS MANUALLY DELETE!!!!!! ####
  input_file <- paste0(rds@project.name,"-condition-DEGs")
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
  
  ## Condition Pathway ####
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
}

# Subcluster to whole ####
Idents(dressSpatial)<- "annotations"
idents <- levels(Lymphocyte)
idents
for (i in idents) {
  select <- WhichCells(Lymphocyte, idents=  i)
  Idents(dressSpatial, cells = select)<- paste("Lymph-", i)
}

DimPlot(dressSpatial,label = T)

dressSpatial$annotations <- Idents(dressSpatial)

# heatmap spatial ####
spatial_metadata <- FetchData(object = seu_obj, vars = c("ident", "imagecol", "imagerow"))

library(ComplexHeatmap)

# Assuming 'proximity_matrix' is an NxN or CellTypeA x CellTypeB matrix
Heatmap(proximity_matrix, 
        name = "Proximity Score", 
        col = colorRampPalette(c("blue", "white", "red"))(50),
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        rect_gp = gpar(col = "white", lwd = 1))

