# Allen Brain Map: MapMyCells
# Generate h5 files under 2G from rds raw counts
library(Seurat)
library(cowplot)
library(reticulate)
library(anndata)
library(rhdf5)
#remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
library(future)
library(hdf5r)
library(readxl)
library(dplyr)
library(purrr)
library(fs)
library(ggplot2)


# Alternative Using zellkonverter ############ 
#BiocManager::install("zellkonverter")
library(zellkonverter)
sce <- as.SingleCellExperiment(rds, assay = c("RNA"))
writeH5AD(sce, paste0(rds@project.name, ".h5ad"), X_name = 'counts')


# Collect Annotations & Map to RDS ####
# Randomly sample n cells
n = 40000
set.seed(99)
selected <- sample(x = colnames(rds), size = n, replace = FALSE)
sub <- subset(rds, cells = selected)
sub@project.name <- paste0(rds@project.name, "-subset", n)
saveRDS(sub, paste0(sub@project.name, ".rds"))
rds <- sub
SaveH5Seurat(rds, filename = paste0(rds@project.name), overwrite = T)
Convert(paste0(rds@project.name, ".h5seurat"),dest = "h5ad", overwrite = T)

# Recursive function to search for MapMyCells excel files and extract cols A & L
mmcExtract <- function(project.name, dir) {
  # initialize list to store all extracted cols
  all_cols = data.frame(cell_id = character(), cluster_name = character(),
                        stringsAsFactors = F)
  
  files <- dir_ls(dir, recurse = T, regexp = "\\.csv$")
  if (length(files) == 0) { stop("No csv files found in directory")}
  print(files)
  print("Starting extract.........")
  
  # Read excel files, ignore first 4 rows of metadata, row 5 is header
  for (file in files){
    data <- read.csv(file, skip = 4, header = T) # skip first 4 rows
    two_cols <- data.frame(cell_id = data[,1], cluster_name = data[,12],
                           stringsAsFactors = F)
    
    # Extract cols A & L
    all_cols = rbind(all_cols, two_cols)
  }
  # Write to csv file
  write.csv(all_cols, paste0(project.name,"-mmcAnno.csv"), row.names = F)
  message(paste("Saved to:", paste0(getwd(), "/",project.name,"-mmcAnno.csv")))
  rm(all_cols)
}
mmcExtract(rds@project.name, paste0(getwd()))

anno <- read.csv(paste0(rds@project.name, "-mmcAnno.csv"), header = T, stringsAsFactors = F)
head(anno)

# Check if all cell IDs in the CSV are present in the Seurat object
# If not, may be due changed cell ID during merge (ogcellid_1 for C001, ogcellid_2 for C002, etc)
# Manually add tag to -mmcAnno.csv

# Check annotations and take off anything with <50? cells
table <- data.frame(table(anno$cluster_name))
table <- table %>% filter(Freq > 50)
table

# Take these small clusters out of anno 
anno <- anno[anno$cluster_name %in% table$Var1,]
table(anno$cluster_name)

all(anno$cell_id %in% Cells(rds))

# Check for mismatches
print(anno$cell_id == Cells(rds))

# If list is different, create new anno with all IDs
MMC <- rep(NA, length(Cells(rds)))
names(MMC) <- Cells(rds)
MMC[anno$cell_id] <-  anno$cluster_name

# Set to rds# Set MMCto rds
rds <- AddMetaData(rds, MMC, col.name = "MMC")
table(rds$MMC)

#PRINT!!
pdf(paste0(rds@project.name,"-UMAPMMC.pdf"), height = 5, width = 10)
DimPlot(rds, reduction = rds@misc$umap, group.by = "MMC")
dev.off()

