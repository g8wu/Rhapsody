# ADTnorm ####
# Tutorial: https://yezhengstat.github.io/ADTnorm/articles/ADTnorm-tutorial.html
BiocManager::install("flowCore")
remotes::install_github("yezhengSTAT/ADTnorm", build_vignettes = F)
library(ADTnorm)

# Get raw ADT counts with cell ID rows, Abseq features cols
cell_x_adt <- t(rds@assays$ADT@counts)
head(cell_x_adt)

# Get desired idents 
cell_x_feature <- data.frame(rds@active.ident)

# Make sure there is a 'sample' column set to an ident
colnames(cell_x_feature) <- 'sample'
head(cell_x_feature)

# Run ADTnorm
cell_x_adt_norm = ADTnorm(
  cell_x_adt = cell_x_adt, 
  cell_x_feature = cell_x_feature, 
  save_outpath = getwd(), 
  study_name = rds@project.name,
  marker_to_process = newAbseq,
  bimodal_marker = NULL,             # default NULL: try different settings to find bimodal peaks
  trimodal_marker = c("CD45"),       # CD4 and CD45RA tend to have 3 peaks
  # setting the CD3 uni-peak of buus_2021_T study to positive peak if only one peak is detected for CD3 marker
  # positive_peak = list(ADT = "CD3", sample = "buus_2021_T"), 
  positive_peak = list(ADT = "CD19", sample = "B Cell"), 
  brewer_palettes = "Dark2",
  save_fig = TRUE,
  target_landmark_location = "fixed",
  shoulder_valley = T,               # Look for "shoulder" as pos peak (technical variation -> no clear separation b/w neg/pos)
  #multi_sample_per_batch = T,        # Omit aligning the one pos peak
  #customize_landmark = T             # Manual adjustment UI 
)

## Put ADTNorm matrix back into original seurat object ####
rds@assays$ADT@data <- t(cell_x_adt_norm)
