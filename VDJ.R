# VDJ ####
# Tutorial: https://www.bioconductor.org/packages/release/bioc/vignettes/scRepertoire/inst/doc/vignette.html
if(!require("scRepertoire")) {
  remotes::install_github("ncborcherding/immApex")
  remotes::install_github("ncborcherding/scRepertoire")
}
n
if(!require("BiocStyle")) {BiocManager::install("BiocStyle")}
n
if(!require("scater")) {BiocManager::install("scater")}
n
library(BiocStyle)
library(scater)
library(Seurat)
library(igraph)
library(scRepertoire)
library(dplyr)

# Read in contig files ####
rds <- readRDS("DressPBMC-BF-ADT-CLR.rds")
S1 <- read.delim("C104-exact-poly10-51x71_VDJ_Dominant_Contigs_AIRR.tsv")
S2 <- read.delim("_1_C105-exact-poly10-51x71_VDJ_Dominant_Contigs_AIRR.tsv")
S3 <- read.delim("C114-expected-poly10-51x71_VDJ_Dominant_Contigs_AIRR.tsv")

contig.list <- list(S1)
contig.list <- loadContigs(input = contig.list, format = "BD")

# Separate TCR and BCR ####
## BCR --------------------------------
airrBCR <- subset(S1, locus %in% c("IGH", "IGK", "IGL"))
contigBCR <- list(airrBCR) %>% loadContigs(format = "BD")
combinedBCR <- combineBCR(contigBCR, samples = "", ID = "")

# strip the 2 underscores added from previous line
combinedBCR[[1]]$barcode <- sub("^_+",  "", combinedBCR[[1]]$barcode)

# Add sample tag metadata
combinedBCR[[1]]$sample <- rds$Sample_Tag[match(combinedBCR[[1]]$barcode, colnames(rds))]
table(combinedBCR[[1]]$sample)

clonalHomeostasis(combinedBCR)
clonalDiversity(combinedBCR)
scRepB <- combineExpression(combinedBCR, rds, cloneCall = "CTstrict")
clonalOverlay(scRepB, reduction = "umap", cut.category = "clonalFrequency") +
  ggtitle("BCR")


## TCR --------------------------------
airrTCR <- subset(S1, locus %in% c("TRA", "TRB", "TRD"))
contigTCR <- list(airrTCR) %>% loadContigs(format = "BD")
combinedTCR <- combineTCR(contigTCR, samples = "", ID = "")

# strip the 2 underscores added from previous line
combinedTCR[[1]]$barcode <- sub("^_+",  "", combinedTCR[[1]]$barcode)

# Add sample tag metadata
combinedTCR[[1]]$sample <- rds$Sample_Tag[match(combinedTCR[[1]]$barcode, colnames(rds))]
table(combinedTCR[[1]]$sample)

clonalHomeostasis(combinedTCR)
clonalDiversity(combinedTCR)
scRepT <- combineExpression(combinedTCR, rds, cloneCall = "CTstrict")
clonalOverlay(scRepT, reduction = "umap", cut.category = "clonalFrequency") +
  ggtitle("TCR")


# QC Plots ####
DimPlot(rds, group.by = "BCR_Paired_Chains")
DimPlot(rds, group.by = "TCR_Paired_Chains")
FeaturePlot(rds, features = "Total_VDJ_Read_Count")

# ClonalCompare Barplot ####
table(combinedBCR[[1]]$sample)
clonalCompare(combinedBCR, 
              top.clones = 100, 
              samples = c("SampleTag10_hs", "SampleTag11_hs", "SampleTag12_hs"), 
              cloneCall="CTgene", 
              graph = "alluvial")

clonalCompare(combinedTCR, 
              top.clones = 10, 
              samples = c("SampleTag10_hs", "SampleTag11_hs", "SampleTag12_hs"), 
              cloneCall="CTstrict", 
              graph = "alluvial")

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
