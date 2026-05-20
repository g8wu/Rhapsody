# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# SET WORK DIRECTORY TO SPECIFIC PROJECT FOLDER
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(Seurat)

sampleSelect <- function(cart, sample, ST, fileType, outDir){
  C1 <- readRDS(paste0(cart, "-", fileType, "_Seurat.rds"))
  
  # If no ST, set entire object as patient
  if(ST != "xx"){
    C1 <- subset(C1, Sample_Tag == paste0("SampleTag", ST))
  }
  
  # Assign orig.ident as patient
  C1@meta.data$orig.ident <- sample
  
  # Subset object and save
  saveRDS(C1, file =  paste0(outDir, "/", sample, "-", fileType, ".rds"))
  print(paste0(cart,"-", ST, " saved as: ", sample, "-", fileType, ".rds to:"))
  print(outDir)
}

#' Reassign ST given threshold pct
#' Input: ...Sample_Tag_ReadsPerCell.csv and respective seurat object
#' @param cutoff pct of total reads needed to be assigned that ST
#' @param fileName name for seurat object and respective .csv file
#' @param STs vector of ints. Expected STs in seurat object
ST.reassign <- function(cutoff, fileName, STs){
  rds <- readRDS(paste0(fileName, "_Seurat.rds"))
  rpc <- read.csv(paste0(fileName, "_Sample_Tag_ReadsPerCell.csv"), skip = 7)
  
  # Adjust ST column index
  STs <- STs + 1
  STs
  # Calculate new STs
  rpc$newST <- apply(rpc[, STs], 1, function(row.counts) {
    total <- sum(row.counts)
    if (total < cutoff) {
      return("Undetermined")
    }
    pcts <- row.counts / total
    
    if (any(pcts >= cutoff)) {
      dom.ST <- names(pcts)[which.max(pcts)]
      return(substr(dom.ST, 1, nchar(dom.ST) - 11))
    } else {
      return("Multiplet")
    }
  })
  
  # Print original ST table
  print("Original ST: 0.75")
  table(rds$Sample_Tag)
  
  # Make seurat object Ident the cell ids
  Idents(rds) <- Cells(rds)
  
  # Reassign STs and save
  newST <- setNames(rpc$newST, rpc[, 1])
  rds <- RenameIdents(rds, newST)
  rds$Sample_Tag <- Idents(rds)
  print(paste("Reassigned ST: ", cutoff))
  table(rds$Sample_Tag)
  saveRDS(rds, paste0(fileName, "_ST.rsgn.", cutoff, ".rds"))
}

# MAIN ##############################

## ST reassign ####
cut <- 0.8

ST.reassign(cut, "C77-exact-poly10-51x71", c(4, 5))
ST.reassign(cut, "C78-exact-poly10-51x71", c(4, 5))
ST.reassign(cut, "C79-exact-poly10-51x71", c(4, 5))
ST.reassign(cut, "C80-exact-poly10-51x71", c(6, 7))
ST.reassign(cut, "C81-exact-poly10-51x71", c(6, 7))
ST.reassign(cut, "C82-exact-poly10-51x71", c(6, 7))
ST.reassign(cut, "C83-exact-poly10-51x71", c(8, 9))
ST.reassign(cut, "C84-exact-poly10-51x71", c(8, 9))

## Separate by ST ####
#mouse "01_mm"
#human "01_hs"
out<- "/mnt/BioAdHoc/Groups/Collaborators/ben.croker/dress/samples/pbmc"

sampleSelect("C86", "DRS03-act-BF", "04_hs", "exact-poly10-51x71", out)


sampleSelect("C85", "DRS02-act", "03_hs", "exact-poly10-51x71", out)
sampleSelect("C93", "DRS02-reco", "01_hs", "exact-poly10-51x71", out)

sampleSelect("C102", "DRS04-act", "01_hs", "exact-poly10-51x71", out)
sampleSelect("C102", "H04", "12_hs", "exact-poly10-51x71", out)
