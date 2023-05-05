setwd('/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2022-02-01_Killifish_scRNAseq_ATLAS_analyses/AUGUR')
options(stringsAsFactors = F)

# load libraries for analysis
library('Seurat')
library(bitops)
library('ggplot2')
library('dplyr')
library(Augur)
library(viridis)

######################################################################
# Load annotated object
load('../Seurat_merged_analysis/2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData')
killi.combined.filt
# An object of class Seurat 
# 38740 features across 81357 samples within 3 assays 
# Active assay: SCT (18370 features, 5000 variable features)
# 2 other assays present: RNA, integrated
# 2 dimensional reductions calculated: pca, umap 

######################################################################
# Augur
augur.killi <-  calculate_auc(as.matrix(killi.combined.filt@assays$SCT@data),
                              killi.combined.filt@meta.data, 
                              cell_type_col = "Annotation_v1", 
                              label_col = "Sex",
                              n_threads = 1)
augur.killi

as.data.frame(augur.killi$AUC)
#                    cell_type       auc
# 1                 Endothelial 0.7916289
# 2                 Hepatocytes 0.6712868
# 3       NK_T_progenitor_cells 0.6462302
# 4  Hepatocytes_Efferocytosing 0.6410488
# 5          Kidney_prox_tubule 0.6050265
# 6                  NK_T_cells 0.5713889
# 7                       HSPCs 0.5705990
# 8                 Fibroblasts 0.5644558
# 9        Kidney_distal_tubule 0.5601115
# 10         B_Cell_Progenitors 0.5518216
# 11                Neutrophils 0.5508258
# 12               Thrombocytes 0.5490495
# 13    Multipotent_progenitors 0.5440665
# 14             Cholangiocytes 0.5378118
# 15        Myeloid_progenitors 0.5278723
# 16     Neutrophil_Progenitors 0.5230707
# 17                    B_cells 0.5210317
# 18    Erythrocyte_Progenitors 0.5185563
# 19                Macrophages 0.5166327
# 20                 Mast_cells 0.5141327
# 21       Lymphoid_progenitors 0.5129592
# 22               Erythrocytes 0.4976096

save(augur.killi, file = paste0(Sys.Date(),"_Augur_Killi_Atlas_object.RData"))

# for some reason, barcode column trips Augur up
killi.combined.filt@meta.data <- killi.combined.filt@meta.data[,colnames(killi.combined.filt@meta.data)  != "barcode"]

pdf(paste0(Sys.Date(),"_Augur_Killi_Atlas_UMAP.pdf"), width = 3, height = 3)
plot_umap(augur.killi,killi.combined.filt, cell_type_col = "Annotation_v1")
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Killi_Atlas_UMAP_Red_Blue.pdf"), width = 3, height = 3)
plot_umap(augur.killi,killi.combined.filt, cell_type_col = "Annotation_v1", palette = colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50))
dev.off()

pdf(paste0(Sys.Date(),"_Augur_Killi_Atlas_Lollipop.pdf"), width = 3, height = 3)
plot_lollipop(augur.killi)
dev.off()
######################################################################



#######################
sink(file = paste(Sys.Date(),"_Killi_Atlas_AUGUR_session_Info.txt", sep =""))
sessionInfo()
sink()

