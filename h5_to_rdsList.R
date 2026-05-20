install.packages("hdf5r")
1
library(Seurat)
library(purrr)

# Get all .h5 files from current directory
wkdir <- getwd()
rds_files <- list.files(path = paste0(wkdir, "/samples/PMID 39814731"), pattern = "h5", full.names = TRUE)
print(rds_files)
list <- list()
list <- map(rds_files, ~ {
  # read h5 file
  matrix <- Read10X_h5(.x)
  # convert to seurat object
  CreateSeuratObject(counts = matrix, min.cells = 3, min.features = 200)
})

print(list)
