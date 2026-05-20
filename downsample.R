# Assuming 'large.obj' is your Seurat object
set.seed(99)  # For reproducibility
Idents(rds) <- rds$orig.ident
wt <- subset(rds, idents = "WT-BRAIN")
nm <- subset(rds, idents = "NOMID-BRAIN")
wtDown <- sample(colnames(wt), size = 49000, replace = FALSE)
nmDown <- sample(colnames(nm), size = 49000, replace = F)

# put together wtDown and nmDown?
downRds <- rds[, c(wtDown, nmDown)]
downRds
table(downRds$orig.ident)

DimPlot(downRds, reduction = "rna.umap")
rds <- downRds
rm(downRds)
rm(nm)
rm(wt)
