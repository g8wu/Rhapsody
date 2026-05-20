# Pseudobulk DEG analysis
library(Seurat)
library(edgeR)
library(dplyr)
library(DESeq2)

# Set group ####
group <-"geno"
#group <- "geno.basicAnno"

Idents(rds) <- group

# Get raw counts
counts <- GetAssayData(rds, assay = "RNA", slot = "counts")
meta <- rds@meta.data
meta$group <- paste(meta$geno, meta$basicAnno, sep = "_")

# Psuedobulk: Aggregate counts/cluster ####
pseudoCounts <- sapply(unique(meta$group), function(g){
  cellIds <- WhichCells(rds, idents = g)
  rowSums(counts[, cellIds, drop=F])
})
pseudoCounts <- do.call(cbind, pseudoCounts)
#pseudoCounts <- as.matrix(pseudoCounts)
colnames(pseudoCounts) <- unique(meta$group)
meta <- data.frame(group = clusters)

# Split 
meta$geno <- ifelse(grepl("WT", meta$group), "WT", "NM")
meta$anno <- sub("_.*", "", meta$group)
rownames(meta) <- colnames(pseudoCounts)

# Create DGE list
dge <- DGEList(counts = pseudoCounts, samples = meta)

# Filter + Norm
keep <- filterByExpr(dge, group = dge$samples$geno)
dge <- dge[keep, , keep.lib.sizes=F]
dge <- calcNormFactors(dge)

# Split by group and model
design <- model.matrix(~0 + clusters, data = meta)
colnames(design) <- clusters

# Dispersion, fit, contrast 
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)

# Choose contrasts
contrast <- makeContrasts(Cluster1-Cluster2, levels = design)
lrt <- glmLRT(fit, contrast = contrast)

topGenes <- topTags(lrt, n = 20)
print(topGenes)