abseq <- rownames(rds[["ADT"]])
abseq

FeatureScatter(rds, feature1 = "Gapdh", feature2 = "adt_CD45")
FeatureScatter(rds, feature1 = "adt_CD11b", feature2 = "adt_CD45")
FeatureScatter(rds, feature1 = "adt_CD11b", feature2 = "adt_CD45")

neuts <- readRDS("C:/Users/gio8w/OneDrive - University of California, San Diego Health/Shared Documents - Rhapsody/Data Analysis/Nomid/Nomid/subcluster/Nomid-RNA0.8-Anno-Neut-RNAClust.rds")
FeatureScatter(neuts, feature1 = "adt_CD11b", feature2 = "adt_CD45")

Idents(neuts) <- neuts$cluster.tissue
FeatureScatter(brainNeut, feature1 = "adt_Ly6G", feature2 = "adt_CD184")
