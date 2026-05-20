# VDJ ####
# Tutorial: https://www.bioconductor.org/packages/release/bioc/vignettes/scRepertoire/inst/doc/vignette.html
if(!require("scRepertoire")) {BiocManager::install("scRepertoire")}
n
if(!require("BiocStyle")) {BiocManager::install("BiocStyle")}
n
if(!require("scater")) {BiocManager::install("scater")}
n
library(BiocStyle)
library(scater)
library(dplyr)
library(Seurat)
library(igraph)
library(scRepertoire) # version 2.x
library(tidyr)
remotes::install_github("corybrunson/ggalluvial")
library(ggalluvial) 
# Downgrade ggplot2
install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz", repos = NULL, type = "source")
library(ggplot2) # version 3.4.4

# .rs.restartR()

# Read in contig files ####
rds <- readRDS("/sbgenomics/project-files/1_Analysis/DressPBMC-BF-WNN.rds")
S1 <- read.delim("/sbgenomics/project-files/C104/exact-poly10-51x71/C104-exact-poly10-51x71_VDJ_Dominant_Contigs_AIRR.tsv")
S2 <- read.delim("/sbgenomics/project-files/C105/exact-poly10-51x71/_1_C105-exact-poly10-51x71_VDJ_Dominant_Contigs_AIRR.tsv")
S3 <- read.delim("/sbgenomics/project-files/C114/exact-poly10-51x71/C114-expected-poly10-51x71_VDJ_Dominant_Contigs_AIRR.tsv")

# TODO: merged rds has cell IDs with extra 1_, 2_, etc. manually add back to the contig data
# "C105-exact-poly10-51x71" suffix _1
# "C104-exact-poly10-51x71" suffix _2
# "DRS03-act-BF"  suffix _3
Idents(rds) <- rds$orig.ident
unique(Idents(rds))
st <- subset(rds, idents = "C105-exact-poly10-51x71")
colnames(st)

S1$cell_id <- paste0(S1$cell_id, "_2")
S2$cell_id <- paste0(S2$cell_id, "_1")
S3$cell_id <- paste0(S3$cell_id, "_3")

contig.list <- list(S1, S2)
# Separate TCR and BCR ####
## BCR --------------------------------
airrBCR.list <- lapply(contig.list, function(x) {
  subset(x, locus %in% c("IGH", "IGK", "IGL"))
})
contig.list <- loadContigs(input = airrBCR.list, format = "BD")
contigBCR <- airrBCR.list %>% loadContigs(format = "BD")
combinedBCR <- combineBCR(contigBCR, 
                          samples = "", 
                          ID = "")

# Add metadata back in
combinedBCR <- lapply(combinedBCR, function(x){
  # Match indces
  indeces <- match(x$barcode, colnames(rds))
  
  # Transfer condition metadata
  x$condition <- rds$condition[indeces]
  x$patient <- rds$patient[indeces]
  x$pat.cond <- rds$pat.cond[indeces]
  
  # Remove NAs 
  x$barcode <- sub("^__",  "", x$barcode)
  x$barcode <- gsub("NA_", "", x$barcode)
  x$CTstrict <- sub("^NA.NA_", "", x$CTstrict)
  x$CTaa <- sub("^NA_", "", x$CTaa)
  
  x
})


# add to seurat object
rds <- combineExpression(combinedBCR, rds, cloneCall = "CTstrict")

pdf(paste0(rds@project.name, "-VDJBQC.pdf"))
clonalHomeostasis(combinedBCR)
clonalDiversity(combinedBCR, group.by = "condition")
scRepB <- combineExpression(combinedBCR, rds, cloneCall = "CTstrict")
clonalOverlay(scRepB, reduction = "rna.umap", cut.category = "clonalFrequency") +
  ggtitle("BCR")
dev.off()

## TCR --------------------------------
contig.list <- list(S1, S2)

airrTCR.list <- lapply(contig.list, function(x) {
  subset(x, locus %in% c("TRA", "TRB", "TRG", "TRD"))
})
contig.list <- loadContigs(input = airrTCR.list, format = "BD")
contigTCR <- airrTCR.list %>% loadContigs(format = "BD")
combinedTCR <- combineTCR(contigTCR, samples = "", ID = "")


# Add metadata back in
combinedTCR <- lapply(combinedTCR, function(x){
  # Match indces
  indeces <- match(x$barcode, colnames(rds))
  
  # Transfer condition metadata
  x$condition <- rds$condition[indeces]
  x$patient <- rds$patient[indeces]
  x$pat.cond <- rds$pat.cond[indeces]
  
  # Remove NAs 
  x$barcode <- sub("^__",  "", x$barcode)
  x$barcode <- gsub("NA_", "", x$barcode)
  x$CTstrict <- sub("^NA.NA_", "", x$CTstrict)
  x$CTaa <- sub("^NA_", "", x$CTaa)
  
  x
})


# add to seurat object
rds <- combineExpression(combinedTCR, rds, cloneCall = "CTstrict")

pdf(paste0(rds@project.name, "-VDJTQC.pdf"))
clonalHomeostasis(combinedTCR)
clonalDiversity(combinedTCR, group.by = "condition")
scRepB <- combineExpression(combinedTCR, rds, cloneCall = "CTstrict")
clonalOverlay(scRepB, reduction = "rna.umap", cut.category = "clonalFrequency") +
  ggtitle("TCR")
dev.off()

# QC Plots ####
pdf(rds@project.name, "-VDJQC.pdf")
DimPlot(rds, group.by = "BCR_Paired_Chains", reduction = rds@misc$umap)
DimPlot(rds, group.by = "TCR_Paired_Chains", reduction = rds@misc$umap)
FeaturePlot(rds, features = "Total_VDJ_Read_Count", reduction = rds@misc$umap)
dev.off()

# Alluvial ####
n <- 10
pdf(paste0(rds@project.name, "-VDJ", n, ".pdf"))
clonalCompare(
  combinedBCR, top.clones = n,
  cloneCall = "CTaa",
  group.by = "condition",
  graph = "alluvial",
  palette = "Dynamic",
  order.by = c("act", "trm", "rec")
) +
  ggtitle(paste("BCR AA top", n))

clonalCompare(
  combinedBCR, top.clones = n,
  cloneCall = "CTaa",
  group.by = "condition",
  graph = "alluvial",
  palette = "Dynamic",
  order.by = c("act", "trm", "rec")
) +
  ggtitle(paste("TCR AA top", n))

dev.off()

## Strict ####
n <- 10
pdf(paste0(rds@project.name, "-VDJstrict", n, ".pdf"))
clonalCompare(
  combinedBCR, top.clones = n,
  cloneCall = "CTaa",
  group.by = "condition",
  graph = "alluvial",
  palette = "Dynamic",
  order.by = c("act", "trm", "rec")
) +
  ggtitle(paste("BCR CTstrict top", n))

clonalCompare(
  combinedBCR, top.clones = n,
  cloneCall = "CTaa",
  group.by = "condition",
  graph = "alluvial",
  palette = "Dynamic",
  order.by = c("act", "trm", "rec")
) +
  ggtitle(paste("TCR CTstrict top ", n))

dev.off()

# # Clone length
# clonalLength(combinedTCR, 
#              cloneCall="aa", 
#              chain = "TRA", 
#              scale = TRUE, 
#              group.by = "sample") 
# 
# # Repertoire metrics ####
# percentAA(combinedTCR, 
#           chain = "TRB", 
#           aa.length = 20)
# 
# positionalEntropy(combinedTCR, 
#                   chain = "TRB", 
#                   aa.length = 20, 
#                   group.by = "sample")
# 
# # The vizGenes() function offers a highly adaptable approach to visualizing the 
# # relative usage of TCR or BCR genes. It acts as a versatile alias for percentGeneUsage(), 
# # allowing for comparisons across different chains, scaling of values, and 
# # selection between bar charts and heatmaps.
# vizGenes(combinedTCR,
#          x.axis = "TRBV",
#          y.axis = NULL, # No specific y-axis variable, will group all samples
#          plot = "barplot", 
#          group.by = "sample") 
# 
# # vizGenes() is particularly useful for examining gene pairings
# vizGenes(combinedTCR,
#          x.axis = "TRBV",
#          y.axis = "TRBJ",
#          plot = "heatmap")
# 
