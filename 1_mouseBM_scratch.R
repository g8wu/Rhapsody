# TODO select RNA umap res 0.3
if(!require("openxlsx", quietly = T)) install.packages("openxlsx")
if(!require("readxl", quietly = T)) install.packages("readxl")
if(!require("EnhancedVolcano", quietly = T)) BiocManager::install("EnhancedVolcano")
library(openxlsx)
library(openxlsx)
library(readxl)
library(patchwork)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(patchwork)
library(EnhancedVolcano)

# Sample naming ####
rds <- readRDS("/sbgenomics/workspace/fatmice/BM/mouseBM-RNAADT.rds")
setwd("/sbgenomics/workspace/fatmice/BM")
rds$geno <- case_when(
  rds$orig.ident == "C111-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag07_mm" ~ "WT",
  rds$orig.ident == "C111-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_mm" ~ "GCSF",
  rds$orig.ident == "C111-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_mm" ~ "IL6",
  rds$orig.ident == "C111-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_mm" ~ "DKO",
  rds$orig.ident == "C112-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag07_mm" ~ "WT",
  rds$orig.ident == "C112-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_mm" ~ "GCSF",
  rds$orig.ident == "C112-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_mm" ~ "IL6",
  rds$orig.ident == "C112-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_mm" ~ "DKO",
  rds$orig.ident == "C113-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag09_mm" ~ "WT",
  rds$orig.ident == "C113-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag10_mm" ~ "GCSF",
  rds$orig.ident == "C113-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag11_mm" ~ "IL6",
  rds$orig.ident == "C113-exact-poly10-51x71" & rds$Sample_Tag == "SampleTag12_mm" ~ "DKO",
  TRUE ~ "Undetermined" # Preserve existing annotations for cells that don't match
)
table(rds$geno)

# set res ####
rds$seurat_clusters <- rds$SCT_snn_res.0.3
Idents(rds) <- "seurat_clusters"

# print Abseq ####
group <- "seurat_clusters"

abseq <- rownames(rds[["ADT"]])
DefaultAssay(rds) <- "ADT"
pdf(paste(rds@project.name, group, "AbseqDot.pdf", sep = "-"),width = 10, height = 6)
DotPlot(rds, features = abseq, col.min = 0, cols = "RdYlBu", group.by = group) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(hjust = 0)) + 
  ggtitle("AbSeq") + 
  #coord_flip() +
  geom_hline(yintercept = seq(n+0.5, length(unique(Idents(rds))) - 0.5, by = n), color = "black")
#geom_hline(yintercept =seq(n-floor(n/2)+ 0.5, length(unique(Idents(rds))) - 0.5, by = n), linetype = "dashed")
#geom_hline(yintercept = custom, color = "black")
dev.off()
DefaultAssay(rds) <- "SCT"

## Featplot #### 
abseq <- rownames(rds@assays$ADT)
DefaultAssay(rds) <-"ADT"
plots <- FeaturePlot(rds, reduction = rds@misc$umap, features = abseq, ncol = 5) & 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank())
DefaultAssay(rds) <-"SCT"

# PRINT
pdf(paste(rds@project.name, "AbsFeats.pdf", sep = "-"), width = length(abseq) / 5 * 4.8, height = length(abseq) / 5 * 5)
print(plots+ plot_annotation(title ="Abseq", theme = theme(plot.title = element_text(size =50))))
dev.off()

degs <- FindAllMarkers(rds, only.pos = F, logfc.threshold = 0.5, assay = "SCT", recorrect_umi = F)
dim(degs)
write.csv(degs, paste0(rds@project.name, "-", group, "-DEGs.csv"))

# Create Excel workbook
wb <- createWorkbook()
for (i in idents){
  print(i)
  subset <- degs[degs$cluster == i,]
  # Sort by decreasing avg LogFC
  subset <- subset[order(subset$avg_log2FC, decreasing = T),]
  addWorksheet(wb, i)
  writeData(wb, sheet = i , subset)
}
saveWorkbook(wb, file=paste0(rds@project.name, "-", group, "-DEGs.xlsx"), overwrite = TRUE)


## Filter DEGs ####
# by LogFC, gate out upregs if pct1 < 0.5
### IF BLANK SHEETS MANUALLY DELETE!!!!!!
input_file <- paste0(rds@project.name, "-", group, "-DEGs")
group <- "basicAnno"
Idents(rds) <- group

# Also dotplot filtered genes
sheet_names <- excel_sheets(paste0(input_file, ".xlsx"))

wb <- createWorkbook()
colors <- colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))
pdf(paste(input_file, "DotTop100.pdf", sep = "-"), width = 20, height = 6)
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
  
  # Print top 100 supervised Dotplot
  print(sheet)
  plot.new()
  # x and y coordinates are from 0 (left/bottom) to 1 (right/top)
  text(x = 0.5, y = 0.5, paste(sheet, "100 Upreg"), cex = 5, font = 2) # cex changes text size, font changes style
  
  Clustered_DotPlot(rds, features = pull(top, gene), flip = T, 
                    colors_use_exp = rev(colors(20)), x_lab_rotate = T, 
                    group.by = group, exp_color_min = 0, 
                    plot_km_elbow = F)
  
  all <- bind_rows(top, bot)
  # Add to workbook
  addWorksheet(wb, sheetName = sheet)
  writeData(wb, sheet = sheet, all)
}
saveWorkbook(wb, file=paste0(input_file, "-filtered.xlsx"), overwrite = TRUE)
dev.off()

## DEG.xlsx -> Volcano ####
#input_file_list <- c("T Cell-GCSF-IL6-DEGs-wilcox", "Monocyte-GCSF-IL6-DEGs-wilcox")

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

