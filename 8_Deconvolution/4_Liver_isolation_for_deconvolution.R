#Subset the completed tissue atlas to isolate liver samples for deconvolution

##########################
#set working directory
setwd("/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/CODE/7_Deconvolution")

#load libraries
library(Seurat)

#generate the liver atlas for deconlvution by subsetting the complete tissue atlas
load("/Volumes/Collaboration_Drive/Killi_Atlas_For_Bryan/Cell_Annotation_Analysis/2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData")

onlyliver <- subset(killi.combined.filt, subset = Tissue %in% "Liver")
save(onlyliver, file = "liver_only_for_deconvolution.Rdata")

#######################
sink(file = paste(Sys.Date(),"_liver_isolation_for_deconvolution_Session_info.txt", sep =""))
sessionInfo()
sink()
