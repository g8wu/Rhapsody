library(harmony)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggpubr)


X10B_subsampled <- readRDS("BAM-C35-Original-Subsampled-250M_Seurat.rds")
X10B_subsampled_51x71 <- readRDS("BAM-C35-Original-Subsampled-250M-51x71-TRIMMED_Seurat.rds")
Aviti <- readRDS("BAM-IndexOriginal_Seurat.rds")
Aviti_reclean_reindex <- readRDS("BAM-ReClean-ReIndex_Seurat.rds")
Aviti_reindex <- readRDS("BAM-ReIndexOnly_Seurat.rds")

X10B_default <- readRDS("C35_Seurat.rds")
X10B_default_cutadapt <- readRDS("C35-cutadapt_Seurat.rds")
X10B_default_noIntron <- readRDS("C35-noIntron_Seurat.rds")
X10B_default_cutadapt_noIntron <- readRDS("C35-cutadapt-noIntron_Seurat.rds")
X10B_exact <- readRDS("C35-exact_Seurat.rds")
X10B_exact_cutadapt <- readRDS("C35-exact-cutadapt_Seurat.rds")
X10B_exact_cutadapt_noIntron <- readRDS("C35-exact-cutadapt-noIntron_Seurat.rds")
X10B_exact_noIntron <- readRDS("C35-exact-noIntron_Seurat.rds")
X10B_exact_poly10_51x71 <- readRDS("C35-exact-poly10-51x71_Seurat.rds")

x <- merge(x = Aviti, y = c(Aviti_reclean_reindex,
                            Aviti_reindex, 
                            X10B_default,
                            X10B_subsampled,
                            X10B_subsampled_51x71,
                            X10B_default_cutadapt, 
                            X10B_default_noIntron, 
                            X10B_default_cutadapt_noIntron, 
                            X10B_exact, 
                            X10B_exact_cutadapt,
                            X10B_exact_cutadapt_noIntron, 
                            X10B_exact_noIntron, 
                            X10B_exact_poly10_51x71))

x@meta.data$orig.ident <- c(rep("Aviti", ncol(Aviti)), 
                            rep("Aviti_reclean_reindex", ncol(Aviti_reclean_reindex)),
                            rep("Aviti_reindex", ncol(Aviti_reindex)), 
                            rep("IlluminaX10B_default", ncol(X10B_default)), 
                            rep("IlluminaX10B_default_subsample", ncol(X10B_subsampled)), 
                            rep("IlluminaX10B_default__subsample_51x71", ncol(X10B_subsampled_51x71)),
                            rep("IlluminaX10B_default_cutadapt", ncol(X10B_default_cutadapt)),
                            rep("IlluminaX10B_default_noIntron", ncol(X10B_default_noIntron)),
                            rep("IlluminaX10B_default_cutadapt_noIntron", ncol(X10B_default_cutadapt_noIntron)),
                            rep("IlluminaX10B_exact", ncol(X10B_exact)), 
                            rep("IlluminaX10B_exact_cutadapt", ncol(X10B_exact_cutadapt)), 
                            rep("IlluminaX10B_exact_cutadapt_noIntron", ncol(X10B_exact_cutadapt_noIntron)), 
                            rep("IlluminaX10B_exact_noIntron", ncol(X10B_exact_noIntron)),
                            rep("IlluminaX25B_exact_cutadapt51x71", ncol(X10B_exact_poly10_51x71)))

rm(X10B_default)
rm(X10B_default_cutadapt)
rm(X10B_default_cutadapt_noIntron)
rm(X10B_default_noIntron)
rm(X10B_exact)
rm(X10B_exact_cutadapt)
rm(X10B_exact_cutadapt_noIntron)
rm(X10B_exact_noIntron)
rm(X10B_exact_poly10_51x71)
rm(X10B_subsampled)
rm(X10B_subsampled_51x71)
rm(Aviti)
rm(Aviti_reclean_reindex)
rm(Aviti_reindex)
gc()

# Set QC thresholds
x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-") 
features = c("nFeature_RNA", "nCount_RNA", "percent.mt") 
print(VlnPlot(x, features = features, pt.size = 0.01, ncol = 3))
x <- subset(x, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA > 150 & nCount_RNA < 30000)
print(VlnPlot(x, features = features, pt.size = 0.01, ncol = 3))

# normalize and identify variable features for each dataset independently
DefaultAssay(x) <- "RNA"
x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
x <- ScaleData(x, verbose = FALSE)
x <- RunPCA(x, features = VariableFeatures(x), npcs = 20, verbose = FALSE)
rm(features)
gc()
# Integrate with harmony
x <- RunHarmony(object = x, "orig.ident")

# Violin plots 
genes <- c("ASIP", "LINC01237", "NUP155", "RGS6", "RORA", "SYPL1")
print(VlnPlot(object = x, features = genes, group.by = "orig.ident", 
              stack = TRUE, flip = TRUE) + NoLegend())

