install.packages("reticulate")
1

remotes::install_github("Moonerss/scrubletR", force = T)
3
library(reticulate)
library(scrubletR)

res <- scrublet_R(seurat_obj = rds)

## use specific threshold
res <- scrublet_R(seurat_obj = rds, threshold = .22)

head(res@meta.data)

scrublet_obj = get_init_scrublet(seurat_obj = samples)

## plot histogram
plot_histogram(scrublet_obj)

## update threshold
scrublet_obj = call_doublets(scrublet_obj, threshold = 0.25)

## plot to check again
plot_histogram(scrublet_obj)

## add info to seurat obj
samples[["doublet_scores"]] <- scrublet_obj$doublet_scores_obs_
samples[["predicted_doublets"]] <- scrublet_obj$predicted_doublets_
