######### Seurat VDJ ########
set.seed(99)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(tidyverse)
library(dplyr)
library(Seurat)
if (!require("scRepertoire", quietly = TRUE)) BiocManager::install("scRepertoire")
1

library(scRepertoire)

rds <- readRDS("Dress0.2.rds")

######## Violin plot features ########
features <- c("Total_VDJ_Read_Count", "Total_VDJ_Molecule_Count")
pdf(paste0(rds@project.name, "-VDJQC.pdf"))
plot_data <- data.frame(x = rds$condition, y = rds$Total_VDJ_Read_Count)
print(ggplot(plot_data, aes(x = x, y = y, fill = x)) +
        geom_violin() + theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90)) +
        labs(title = "Total_VDJ_Read_Count", x = "Condition", y = "Total_VDJ_Read_Count"))

######## scRepertoire ########

# Load necessary libraries
library(Seurat)
library(tidyverse)
library(scRepertoire)

# Assuming `rds` contains BD Rhapsody data
# Load V(D)J data (modify path if necessary)
vdj_data <- read.csv("path/to/vdj_data.csv")  # Adjust the file path

# Integrate V(D)J data with the Seurat object
rds <- combineExpression(vdj_data, rds, cloneCall = "gene")

# Explore V(D)J clonotype distribution
clonal_distribution <- clonalProportion(rds, cloneCall = "gene")
print(clonal_distribution)

# Visualize clonotype diversity
clonalViz <- clonalDiversity(rds, cloneCall = "gene", exportTable = TRUE)
print(clonalViz)

# Compare clonotype between clusters
clonal_cluster_plot <- clonalNetwork(rds, cloneCall = "gene", proportion = TRUE)
print(clonal_cluster_plot)

# Save processed Seurat object
saveRDS(rds, file = "processed_rds.rds")

# Output summary
print("V(D)J analysis completed successfully!")