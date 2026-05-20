# QC ####
# Pick and merge all your seurat objects into one, QC threshold
set.seed(99)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(Seurat)

# Project Settings #####
project <- "dressWB"
fileType <- ".rds"
# wkdir <- setwd("/mnt/BioAdHoc/Groups/Collaborators/ben.croker/fatmice/BM")
wkdir  <- getwd()
print("_______________________________________________")

# For laptops/PCs override default variable size limit
#n <- 50     # Let variables reach up to n GB
#options(future.globals.maxSize= n * 1e9)  # x * 1e9 = x GB

# Read in files ####
rds_files <- list.files(path = wkdir, pattern = fileType, full.names = TRUE)
print(rds_files)
cat(rds_files, file = "samples.txt", sep = "\n")
list <- list()
list <- map(rds_files, readRDS)
print(rownames(list[[1]]@assays$ADT))

########### Rename Abseq ###########
# NOMID
# abseq <- c("CD105", "CD115", "CD117", "CD11a", "CD11b", "CD11c", "CD162", "CD16/32",
#            "CD184", "CD19", "CD335", "CD41", "CD45", "CD48", "CD62L", "CD71",
#            "CXCR2", "Clec7a", "CD150", "F4/80", "Ly6A/E", "Ly6G", "NK1.1", "SiglecF",
#            "TCRb", "TER119")
# DRESS
abseq <- c("CD101", "CD10", "CD11b", "CD11c", "CD123", "CD14", "CD15", "CD162-SELPLG",
           "CD16-FCGR3A", "CD183-CXCR3", "CD184-CXCR4", "CD193-CCR3", "CD194-CCR4",
           "CD19", "CD32-FGCR2A", "CD33", "CD34", "CD3", "CD41-ITGA2B", "CD44",
           "C56-NCAM16.2", "CD62L", "CD63", "CD64-FCGR1A", "CD86", "CD95-FAS", "CXCR2",
           "Siglec8", "CD66b", "FCER1A", "HLA-DR-CD74")

# HFD MOUSE
# abseq <- c("CD105", "CD115-Csf1r", "CD117-Kit", "CD11a-Itgal", "CD11b-ITGAM", 
#            "CD11c-Itgax", "CD162-Selplg","CD16-CD32", "CD184-Cxcr4", "CD19", "CD335-Ncr1", 
#            "CD41-Itga2b", "CD45-F11-Ptprc", "CD48", "CD62L-Sell", "CD71-Tfrc", 
#            "CXCR2", "CD150", "F4-80-Adgre1", "Ly6A-Ly6E", "Ly6G", "NK-1.1-Klrb1b-cKlrb1c",
#            "SiglecF", "Tcrb", "TER119-Ly76")

# DUPI
# abseq <- c("C101", "CD10", "CD11b", "CD11c", "CD123-IL3RA", "CD14", "CD15-FUT4", 
#            "CD162-SELPLG", "C16", "CD183-CXCR3", "CD184-CXCR4", "CD193-CCR3", 
#            "CD194-CCR4", "CD19", "CD32", "CD33", "CD34", "CD3", "CD41-ITGA2B", 
#            "CD44", "CD56", "CD62L", "CD63", "CD64", "CD86", "CD95", "CXCR2", 
#            "Siglec8", "CD66b", "FCER1A", "HLA-DR-CD74")


# Abseq Match Check ####
abseqCheck <- sapply(list, function(x){
  all(rownames(list[[1]]@assays$ADT) %in% rownames(x@assays$ADT))
})
if(!all(abseqCheck)){
  print("!!!! ERROR: ABSEQ MISMATCH IN RDS LIST !!!!")
  print(which(!abseqCheck))
} else {print("Abseq match check: PASS")}


# Merge (no batch effect correct) #### 
rds <- merge(x = list[[1]], y = list[-1], project = project)
write.csv(table(rds$orig.ident), paste0(project, "_cartCounts.csv"))

# QC First Look ####
# Calc mito pct and Rename ADT ####
features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "ADT_total")

# Replace Abseq names
temp <- rds@assays$ADT@counts
rownames(temp) <- abseq
rds[["ADT"]] <- CreateAssayObject(counts = temp)
rownames(rds@assays$ADT)
rm(temp)

# Calc total abseq reads
rds$ADT_total <- colSums(rds@assays$ADT@data)

# Calc percent mitochondrial reads
rds[["percent.mt"]] <- PercentageFeatureSet(rds, pattern = "(?i)^mt-")


# PRINT each feature as one Violin plot
rds$cart.ST <- paste(rds$orig.ident, rds$Sample_Tag)
pdf(paste(project, "VlnQC.pdf", sep = "-"), width = 20, height = 7)
print(VlnPlot(rds, features = features, pt.size = 0, ncol = 4, group.by = "orig.ident") + NoLegend())
print(VlnPlot(rds, features = features, pt.size = 0, ncol = 4, group.by = "cart.ST") + NoLegend())
dev.off()

# Apply QC thresholds #####
nFeatLower <- 50
nFeatUpper <- 2500
nCountLower <- 50
nCountUpper <- 10000
mtPct <- 25

cat("nFeatLower <- 50
nFeatUpper <- 2500
nCountLower <- 50
nCountUpper <- 10000
mtPct <- 25", file = "QCthresholds.txt")

rdsPostQC <- subset(rds, subset = nFeature_RNA > nFeatLower & nFeature_RNA < 
                      nFeatUpper & nCount_RNA > nCountLower & nCount_RNA < 
                      nCountUpper & percent.mt < mtPct)

# PRINT
pdf(paste(project, "VlnPostQC.pdf", sep = "-"), width = 15, height = 10)
print(VlnPlot(rdsPostQC, features = features, pt.size = 0, ncol = 4, group.by = "orig.ident") + NoLegend())
print(VlnPlot(rdsPostQC, features = features, pt.size = 0, ncol = 4, group.by = "cart.ST") + NoLegend())
dev.off()

# SAVE POSTQC ####
saveRDS(rdsPostQC, paste0(rds@project.name, ".rds"))
