# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
# ## SET WORK DIRECTORY TO SPECIFIC PROJECT FOLDER
# ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
#library(Seurat)

cart1 <-"C67"
sample <- "C57.6-LUNG-INFECTED"
preST <- "03_mm"

fileType <- "-mouseCocci-exact50k-p10-51x71"
##############################
# Read in object
skibidi <- paste0(getwd())
C1 <- readRDS(paste0(cart1, fileType, "_Seurat.rds"))

# If no ST, set entire object as patient
if(preST != "xx"){
  C1 <- subset(C1, Sample_Tag == paste0("SampleTag", preST))
}

# Assign orig.ident as patient
C1@meta.data$orig.ident <- sample

# Error log output
line1 <- paste0(Sys.time(), "\n", cart1, fileType, "-ST", preST, " >> pre", sample, "\n")

# Subset object and save
tryCatch({
  saveRDS(C1, file =  paste0(skibidi,  "/", sample, fileType, ".rds"))
  write(line1, file = "LOG.txt", append = TRUE)
  print(paste0(sample, fileType, ".rds", " done!"))
}, error = function(e){
  write(paste0("Error occured: ",line1, e$message))
  print("Error occured: ", sample )
})
