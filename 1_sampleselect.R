library(Seurat)
wkdir <- getwd()

pjt <- "cocci-mouseOnly"
ptn <- "C57.1"
cart1 <-"C60"
preST <- "01"
# cart2 <- "C51"
# postST <- "09"

fileType <- "mouseOnly-exact-poly10-51x71"
##############################
C1 <-  readRDS(paste0(wkdir, "/", cart1,"-", fileType, "_Seurat.rds"))
if(preST != "xx"){
  C1 <- subset(C1, Sample_Tag == paste0("SampleTag", preST,"_hs"))
}
C1@meta.data$orig.ident <- paste0("remiss", ptn)
line1 <- paste0(Sys.time(), "\n", cart1, "-", fileType, "-ST", preST, " >> pre", ptn, "\n")
tryCatch({
  saveRDS(C1, file =  paste0(wkdir,  "/", pjt, "/remiss", ptn, "-",fileType, ".rds"))
  write(line1, file = "sampleSelectLog.txt", append = TRUE)
  print("remiss done!")
}, error = function(e){
  write(paste0("Error occured: ",line1, e$message))
  print("Error occured: pre", ptn )
})

#############
# C2 <-  readRDS(paste0(wkdir, "/", cart2, "-", fileType, "_Seurat.rds"))
# if(postST != "xx"){
#   C2 <- subset(C2, Sample_Tag == paste0("SampleTag", postST,"_hs"))
# }
# C2@meta.data$orig.ident <- paste0("post", ptn)
# line2 <- paste0(Sys.time(), "\n", cart2, "-", fileType, "-ST", postST, " >> post", ptn, "\n")
# tryCatch({
#   saveRDS(C2, file =  paste0(wkdir, "/", pjt, "/post", ptn, "-",fileType, ".rds"))
#   write(line2, file = "sampleSelectLog.txt", append = TRUE)
#   print("Post done!")
# }, error = function(e){
#   write(paste0("Error occured: ", line2, e$message))
#   print("Error occured: post", ptn )
# })
