# DEG
if(!require("openxlsx", quietly = T)) install.packages("openxlsx")
library(openxlsx)
library(readxl)
library(Seurat)
library(harmony)
library(DESeq2)
library(dplyr)
library(harmony)
library(tidyverse)
library(SingleCellExperiment)
library(sctransform)
library(ggplot2)
library(purrr)
set.seed(99)

setwd("/mnt/BioAdHoc/Groups/Collaborators/ben.croker/dress")
rds <- PrepSCTFindMarkers(rds)

resRange <- seq(0.8, 0.7, by = 0.05)
algKey <- c("Louvain", "Refined Louvain", "SLM", "Leiden")
alg <- 4
test <- "bimod"
types <- rownames(table(Idents(rds)))
types
group <- "strain"
groups <- list(c("NM", "WT"))

# RUN SUBSETTING & DEG
for (cluster in types){
  sub <- subset(rds, idents = cluster)
  list  <- SplitObject(sub, split.by = "orig.ident")

  # RNA Normalization
  list <- lapply(list, SCTransform, vars.to.regress = "percent.mt", verbose = F)
  
  # RNA Integration
  features <- SelectIntegrationFeatures(list, nfeatures = 3000, verbose = F)
  list <- PrepSCTIntegration(list, anchor.features = features)
  anchors <- FindIntegrationAnchors(list, normalization.method = "SCT", 
                                    anchor.features = features, dims = 1:30, verbose = F)
  
  combinedRNA <- IntegrateData(anchors, normalization.method = "SCT", k.weight = 30)
  combinedRNA <- ScaleData(combinedRNA)
  combinedRNA <- RunPCA(combinedRNA)

  # Batch effect correction
  combinedRNA <- RunHarmony(combinedRNA, group.by.vars = "cart",
                            reduction.use = "pca", reduction.save = "harmony.rna")

  combinedRNA <- RunUMAP(combinedRNA, reduction = "pca", reduction.name = "rna.umap",
                         reduction.key = 'rnaUMAP_', dims = 1:20)

  # ADT Norm
  abseq <- rownames(list[[1]]@assays$ADT)
  list <- lapply(list, function(i) {
    DefaultAssay(i) <- 'ADT'
    i <- NormalizeData(i, normalization.method = 'CLR', margin = 2)
    i <- ScaleData(i, features = abseq)
    return(i)
  })

  # ADT Integration
  assays <- rep(c("ADT"), times = length(list))
  anchors <- FindIntegrationAnchors(list, anchor.features = abseq, assay = assays)
  combinedADT <- IntegrateData(anchors, dims = 1:20, k.weight = 33)

  DefaultAssay(combinedADT) <- "ADT"
  # Only features are from abseq panel
  combinedADT <- FindVariableFeatures(combinedADT, features = rownames(combinedADT[["ADT"]]))
  combinedADT <- ScaleData(combinedADT, verbose = FALSE)
  combinedADT <- RunPCA(combinedADT, reduction.name = 'pca', approx = FALSE)
  # Batch effect correction
  combinedADT <- RunHarmony(combinedADT, group.by.vars = "orig.ident",
                            reduction = "pca", reduction.save = "harmony.adt",
                            reduction.key = 'adtUMAP_', assay.use = "ADT")
  combinedADT <- RunUMAP(combinedADT, reduction = "pca", reduction.name = "adt.umap",
                         reduction.key = 'adtUMAP_', dims = 1:20)

  combinedRNA[["ADT"]] <- combinedADT[["ADT"]]
  combinedRNA[["apca"]] <- combinedADT[["pca"]]
  combinedRNA[["harmony.adt"]] <- combinedADT[["harmony.adt"]]
  combinedRNA[["adt.umap"]] <- combinedADT[["adt.umap"]]

  sub <- combinedRNA
  sub@project.name <- paste0(rds@project.name, "-", cluster)

  sub <- FindVariableFeatures(sub, assay = "RNA", selection.method = "vst", nfeatures =3000)
  sub <- FindNeighbors(sub, dims = 1:30, graph.name = "RNA_nn")
  sub <- RunHarmony(sub, group.by.vars = "orig.ident",
                    reduction.use = "pca", reduction.save = "harmony.rna")
  sub <- RunUMAP(sub, dims = 1:30, reduction.name = "rna.umap")
  sub@misc$subclustUmap <- "rna.umap"
  
  sub@misc$subclustUmap <- "wnn.umap"
  sub@misc$subclustAlg <- algKey[alg]

  # Multi-res Clustering
  pdf(paste0(sub@project.name, "-RNAWNNUMAP.pdf"))
  sub <- FindClusters(sub, algorithm = alg, resolution = res, graph.name = "RNA_nn")
  Idents(sub) <- sub[[paste0("RNA_nn_res.", res)]][,1]
  print(DimPlot(sub, reduction = sub@misc$subclustUmap, label = TRUE, raster = F) +
          labs(title = paste0("RNA WNN UMAP ", sub@misc$clustAlg,": ", res)))
  dev.off()
  
  ## RNA & ADT WNN ############## 
  sub <- FindVariableFeatures(sub, assay = "integrated", selection.method = "vst", nfeatures = 3000)
  sub <- FindMultiModalNeighbors(sub, reduction.list = list("pca", "apca"),
                                 dims.list = list(1:20, 1:20), modality.weight.name = "integrated.weight")
  sub <- RunUMAP(sub, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  # Save metadata 
  #sub@project.name <- paste0(project, -"WNN")
  # sub@misc$clustAlg <- algKey[alg]
  # sub@misc$umap <- "wnn.umap"
  
  sub <- FindClusters(sub, algorithm = alg, resolution = res, graph.name = "wsnn")
  Idents(sub) <- sub[[paste0("wsnn_res.", res)]][,1]
  
  pdf(paste0(sub@project.name, "-WNNUMAP.pdf"))
  print(DimPlot(sub, reduction= sub@misc$subclustUmap, label = T, raster = F) + 
          labs(title = paste0(rds@project.name, "ADT&RNA UMAP ", algKey[alg],": ", res)))

  features <-c("nFeature_RNA", "nCount_RNA", "percent.mt")
  print(VlnPlot(sub, features = features, pt.size = 0, ncol = length(features)))
  for (feat in features){
    print(FeaturePlot(sub, reduction = sub@misc$subclustUmap, features = feat))
  }
  
  # Print by sample
  Idents(sub) <- sub$tissue.geno
  idents <- levels(Idents(sub))
  plots <- for(i in idents){
    cells <- WhichCells(sub, idents = i)
    print(DimPlot(sub, reduction =  sub@misc$subclustUmap, cells.highlight = cells, raster=F) +
            labs(title = i) + NoLegend())
  }
  dev.off()
  
  saveRDS(sub, file= paste0(sub@project.name, "-RNAClust.rds"))
  print(paste0(Sys.time(), " -> RNA ONLY multi res WNN done!"))

  # DEG
  sub <- PrepSCTFindMarkers(sub)
  Idents(sub) <- sub[[group]][,1]
  #table(Idents(sub))
  wb <- createWorkbook()
  for (pair in groups) {
    print(paste(pair[[1]][1], "vs", pair[[2]][1]))
    degs <- FindMarkers(sub, test.use = test, only.pos = F, ident.1 = pair[[1]][1], ident.2 = pair[[2]][1],
                        min.pct = 0.1, logfc.threshold = 0.15, assay = "SCT")
    degs$genes <- rownames(degs)
    degs <- degs[order(degs$avg_log2FC, decreasing = T),]
    addWorksheet(wb, sheet = paste0(pair[[1]][1],".", pair[[2]][1]))
    writeData(wb, sheet = paste0(pair[[1]][1],".", pair[[2]][1]), degs)
  }
  saveWorkbook(wb, file=paste0(sub@project.name, "-", "-DEGs-", test, ".xlsx"), overwrite = TRUE)
}

## DEG w/ already subsetted ####
test<- "bimod"
for (cluster in types){
  sub <- readRDS(paste0(cluster, "-RNAClust.rds"))
  sub <- PrepSCTFindMarkers(sub)
  Idents(sub) <- sub$geno
  wb <- createWorkbook()
  for (pair in groups) {
    print(paste(pair[[1]][1], "vs", pair[[2]][1]))
    degs <- FindMarkers(sub, test.use = test, only.pos = F, ident.1 = pair[[1]][1], ident.2 = pair[[2]][1],
                        min.pct = 0.1, logfc.threshold = 0.15, assay = "SCT")
    degs$genes <- rownames(degs)
    degs <- degs[order(degs$avg_log2FC, decreasing = T),]
    addWorksheet(wb, sheet = paste0(pair[[1]][1],".", pair[[2]][1]))
    writeData(wb, sheet = paste0(pair[[1]][1],".", pair[[2]][1]), degs)
  }
  saveWorkbook(wb, file=paste0(cluster, "-", pair[[1]][1],"-", pair[[2]][1], "-DEGs-", test, ".xlsx"), overwrite = TRUE)
}


## Filter and Volcano ####
# by LogFC, gate out upregs if pct1 < 0.5
input_file_list <- c("Nomid-RNA0.8-Anno-Neuts-cluster.tissue-DEGs")

for(input_file in input_file_list) {
  sheet_names <- excel_sheets(paste0(input_file, ".xlsx"))
  
  wb <- createWorkbook()
  for (sheet in sheet_names) {
    # Read the sheet
    data <- read_excel(paste0(input_file, ".xlsx"), sheet = sheet)
    
    # Apply filtering logic
    filtered <- data %>%
      filter(p_val_adj < 0.05 & (
        (avg_log2FC > 0 & pct.1 >= 0.5) |
          (avg_log2FC < 0 & pct.2 >= 0.5)
      )) %>%
      { rbind(head(., 100), tail(., 100)) }
    
    # Add to workbook
    addWorksheet(wb, sheetName = sheet)
    writeData(wb, sheet = sheet, filtered)
  }
  saveWorkbook(wb, file=paste0(input_file, "-filtered.xlsx"), overwrite = TRUE)
}

## DEG.xlsx -> Volcano ####
#input_file_list <- c("T Cell-GCSF-IL6-DEGs-wilcox", "Monocyte-GCSF-IL6-DEGs-wilcox", "Neutrophil-GCSF-IL6-DEGs-wilcox")

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