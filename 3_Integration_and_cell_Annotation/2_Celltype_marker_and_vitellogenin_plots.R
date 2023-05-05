#Get the median unique genes per cell type per sex and generate pie charts of cell type proportions by sex and tissue

##########################
#set working directory
setwd("/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/CODE/2_Integration_and_cell_Annotation")

#load libraries
library(Seurat)
library(ggpubr)

# get UMAP
load("/Volumes/Collaboration_Drive/Killi_Atlas_For_Bryan/Cell_Annotation_Analysis/2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData")

#1. Plot Vitelogenin-related plots

#Featureplot
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsFeature_plot_Singlets_vitellogenin.pdf", sep = "_"), height = 15, width = 13)
FeaturePlot(killi.combined.filt, features = "LOC107392619", cols = c("cornsilk", "red"), order = T, raster = T)
dev.off()

#Violin plot of cell types by sex
pdf(paste(Sys.Date(),"vitellogenin_celltype_violin.pdf", sep = "_"))
VlnPlot(killi.combined.filt, features = "LOC107392619", group.by = "Annotation_v1",  cols = c("Deeppink", "Deepskyblue"),  pt.size = 0, assay = 'RNA', split.by = "Sex")
dev.off()

#Violin plot of tissue by sex
pdf(paste(Sys.Date(),"vitellogenin_tissue_violin.pdf", sep = "_"))
VlnPlot(killi.combined.filt, features = "LOC107392619",  group.by = c("Tissue", "Annotation_v1"), cols = c("Deeppink", "Deepskyblue"),  pt.size = 0, assay = 'RNA', split.by = "Sex")
dev.off()

################

#2. Plot select markers as featureplots

#set parameters 
colors = c("cornsilk", "red")

#organize markers
bcells <- c("pax5","LOC107387651","cd79a")
macrophages <- c("csf1r","marco","LOC107394065")
neutrophils <- c("LOC107372911","LOC107386182","nccrp1")
hepatocytes <- c("LOC107383970","fgg","LOC107378176")

#plot as concatendated featureplots
bcell_plot <- FeaturePlot(killi.combined.filt, features = bcells, cols = colors, order = T, raster = T )
macrophage_plot <- FeaturePlot(killi.combined.filt, features = macrophages, cols = colors, order = T, raster = T )
neutrophil_plot <- FeaturePlot(killi.combined.filt, features = neutrophils, cols = colors, order = T, raster = T )
hepatocyte_plot <- FeaturePlot(killi.combined.filt, features = hepatocytes, cols = colors, order = T, raster = T )

#print all plots together
pdf(paste(Sys.Date(),"test_ftplot.pdf", sep = "_"), height = 30, width = 9)
ggarrange(bcell_plot, macrophage_plot, neutrophil_plot, hepatocyte_plot, nrow = 4)
dev.off()

#######################
sink(file = paste(Sys.Date(),"_Celltype_marker_and_vitellogenin_plots.txt", sep =""))
sessionInfo()
sink()
