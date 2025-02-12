library(dplyr)
library(Seurat)
library(patchwork)

setwd("~/Research/Rhapsody")

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

file <- "pbmc3k_filtered_gene_bc_matrices.tar.gz"


