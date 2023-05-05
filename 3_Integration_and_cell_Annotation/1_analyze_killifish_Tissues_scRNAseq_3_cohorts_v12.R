setwd('/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2022-02-01_Killifish_scRNAseq_ATLAS_analyses/Seurat_merged_analysis')
options(stringsAsFactors = F)

# Single Cell and general use packages
library('Seurat')    # 
library(bitops)      # 
library(sctransform) # 
library(pheatmap)    # 
library(clustree)    # 
library(scales)      # 
library(dplyr)       # 
library(Vennerable)  # 
library(beeswarm)    # 
library(scales)      #

##########  Cell identity annotation packages ##########  
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
library(scSorter)

# # https://github.com/pcahan1/singleCellNet
library(singleCellNet)

#https://github.com/Irrationone/cellassign
# devtools::install_version("tensorflow", version = "2.4.0")
# tensorflow::install_tensorflow(extra_packages='tensorflow-probability', version = "2.1.0")
# tensorflow::tf_config()
# # Loaded Tensorflow version 2.1.0
# # TensorFlow v2.1.0 (~/.virtualenvs/r-reticulate/lib/python3.7/site-packages/tensorflow)
# # Python v3.7 (~/.virtualenvs/r-reticulate/bin/python)

# devtools::install_github("Irrationone/cellassign")
library(cellassign)              #
library(tensorflow)              #
library(SingleCellExperiment)    #
########################################################


########################################
# 2022-09-09
# Restart with data cleaned with SoupX (remove ambient RNA) and with intersection of doubletFinder and scds

# 2022-09-19
# computer had crashed while at FASEB conference. Resume CellAssign using saved intermediate objects

# 2022-10-04
# Based on extreme "mixing" between cell types by cellassign, the intersection rather than union keeps too many doublets.
# Better to be conservative in this case
# Rerun with SoupX but with union
# Will try scSorter instead of CelllAssign since it may outperform
### see: https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html

# 2022-10-10
# resume run at scSorter step

# 2022-10-11
# Also include independent try of Cross-Species SingleCellNet (MWS reference)
# And try scType
# https://github.com/IanevskiAleksandr/sc-type ; https://www.nature.com/articles/s41467-022-28803-w

# 2022-10-11
# add batch integration (with Seurat RPCA)
# Seurat bug with "split.by" made it seem like there was good mixing
# Actually, batches separate quite a bit

# 2023-01-31
# Bryan did manual curation of unannotated or ambiguous clusters
# Integrate and add to object for analysis

# 2023-02-08
# manually curate annotations and identify markers to transfer

#####################################################################################################################
#### 1. Load Cleaned up Seurat Objects

# load up cleaned up singlets
load("../Clean_Seurat_Singlets/2022-09-25_Killifish_Tissues_Cohort1_Seurat_object_SINGLETS.RData")
load("../Clean_Seurat_Singlets/2022-10-04_Killifish_Tissues_Cohort2_Seurat_object_SINGLETS.RData")
load("../Clean_Seurat_Singlets/2022-10-04_Killifish_Tissues_cohort3_Seurat_object_SINGLETS.RData")

killi.singlets.c1
# An object of class Seurat 
# 36066 features across 41195 samples within 2 assays 
# Active assay: SCT (18033 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

killi.singlets.c2
# An object of class Seurat 
# 34492 features across 20870 samples within 2 assays 
# Active assay: SCT (17246 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

killi.singlets.c3
# An object of class Seurat 
# 31560 features across 19367 samples within 2 assays 
# Active assay: SCT (15780 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

# Merge Seurat objects
killi.singlets.merged <- merge(killi.singlets.c1,
                               y =  c(killi.singlets.c2  ,
                                      killi.singlets.c3  ),
                               project = "10x_killi_tissues")

killi.singlets.merged
# An object of class Seurat 
# 36740 features across 81432 samples within 2 assays 
# Active assay: SCT (18370 features, 0 variable features)
# 1 other assay present: RNA

# bring RNA as main assay again
DefaultAssay(killi.singlets.merged) <- "RNA"

table(killi.singlets.merged@meta.data$Batch)
# Cohort_1 Cohort_2 Cohort_3
#      41195    20870    19367 


# need to clean up and rerun SCT on merged data
# https://github.com/satijalab/seurat/issues/2662
killi.singlets.merged[['SCT']] <- NULL
# An object of class Seurat 
# 18370 features across 81432 samples within 1 assay 
# Active assay: RNA (18370 features, 0 variable features)

rm(killi.singlets.c1,killi.singlets.c2,killi.singlets.c3)

################################################################################################################################################################
#### 2. QC on mitochondrial reads
# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsviolinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 15)
VlnPlot(object = killi.singlets.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "Batch")
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(killi.singlets.merged, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "Batch")
plot2 <- FeatureScatter(killi.singlets.merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Batch")

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsQC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2))
dev.off()

################################################################################################################################################################
#### 3. Store info on biological origin of sample

# remove irrelevant columns (previous SCT on single object, Doublet Finder column and old clustering info)
sct.cols <- grep("SCT",colnames(killi.singlets.merged@meta.data))
rm.cols  <- c(grep("seurat_clusters",colnames(killi.singlets.merged@meta.data)),
              grep("DoubletFinder",colnames(killi.singlets.merged@meta.data))  ,
              grep("scds_hybrid",colnames(killi.singlets.merged@meta.data)))
killi.singlets.merged@meta.data <- killi.singlets.merged@meta.data[,-c(sct.cols,rm.cols)]


killi.singlets.merged@meta.data$Group <- paste(killi.singlets.merged@meta.data$Tissue,killi.singlets.merged@meta.data$Sex, sep = "_")

head(killi.singlets.merged@meta.data)
#                                            orig.ident nCount_RNA nFeature_RNA            Condition Sex Tissue    Batch percent.mito      S.Score    G2M.Score
# Female_Blood_cohort1_AAACCCAAGAAGAACG-1 SeuratProject       2874          487 Female_Blood_cohort1   F  Blood Cohort_1     0.000000  0.042706150  0.016107114
# Female_Blood_cohort1_AAACCCAAGAAGTGTT-1 SeuratProject       1175          351 Female_Blood_cohort1   F  Blood Cohort_1     0.000000 -0.025396066 -0.008228870
# Female_Blood_cohort1_AAACCCAAGCCATTCA-1 SeuratProject       9793         2343 Female_Blood_cohort1   F  Blood Cohort_1     3.696518 -0.034265817 -0.095839758
# Female_Blood_cohort1_AAACCCAGTACCGGCT-1 SeuratProject       6375          430 Female_Blood_cohort1   F  Blood Cohort_1     0.000000  0.004213373  0.034397162
# Female_Blood_cohort1_AAACCCAGTGATGTAA-1 SeuratProject        772          468 Female_Blood_cohort1   F  Blood Cohort_1     0.000000 -0.006531250 -0.122278665
# Female_Blood_cohort1_AAACGAAAGCAGCCTC-1 SeuratProject      24771          357 Female_Blood_cohort1   F  Blood Cohort_1     0.000000  0.005052965 -0.001057473
#                                         Phase     old.ident DoubletCall   Group
# Female_Blood_cohort1_AAACCCAAGAAGAACG-1     S SeuratProject     Singlet Blood_F
# Female_Blood_cohort1_AAACCCAAGAAGTGTT-1    G1 SeuratProject     Singlet Blood_F
# Female_Blood_cohort1_AAACCCAAGCCATTCA-1    G1 SeuratProject     Singlet Blood_F
# Female_Blood_cohort1_AAACCCAGTACCGGCT-1   G2M SeuratProject     Singlet Blood_F
# Female_Blood_cohort1_AAACCCAGTGATGTAA-1    G1 SeuratProject     Singlet Blood_F
# Female_Blood_cohort1_AAACGAAAGCAGCCTC-1     S SeuratProject     Singlet Blood_F


### Cell cycle frequency
freq.cycle <- prop.table(x=table(killi.singlets.merged@meta.data$Phase, killi.singlets.merged@meta.data$Condition ), margin =2)

pdf(paste0(Sys.Date(), "_killi_atlas_cell_cycle_phases_barplot.pdf"), height = 5, width= 10)
barplot( freq.cycle , las = 2)
dev.off()

################################################################################################################################################################
#### 4. Scaling the data and removing unwanted sources of variation

# global-scaling normalization method 
# SCT replaces Normalize/FindVariableFeatures/ScaleData
killi.atlas.clean     <- SCTransform(object = killi.singlets.merged, vars.to.regress = c("nFeature_RNA", "percent.mito", "Batch"), variable.features.n = 5000)
killi.atlas.clean
# An object of class Seurat 
# 36740 features across 81432 samples within 2 assays 
# Active assay: SCT (18370 features, 5000 variable features)
# 1 other assay present: RNA

table(killi.atlas.clean@meta.data$Batch, killi.atlas.clean@meta.data$Group)
#           Blood_F Blood_M Kidney_F Kidney_M Liver_F Liver_M Spleen_F Spleen_M
# Cohort_1    4434    4525     6770     5064    3852    6983     4584     4983
# Cohort_2    2850    3470     4922        0    2797    2713        0     4118
# Cohort_3       0       0     4125     3864    2523    2512     2866     3477

save(killi.atlas.clean, file = paste0(Sys.Date(),"_Killi_Fish_3cohorts_Seurat_object_SINGLETS_postSCT.RData"))


###############################################################################################################################################################
##### 5. Cluster the cells
# load object
load('2022-10-04_Killi_Fish_3cohorts_Seurat_object_SINGLETS_postSCT.RData')

# Run dimensionality reduction PCA
killi.atlas.clean <- RunPCA(killi.atlas.clean, npcs = 50)

# Determine the ‘dimensionality’ of the dataset
# Approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
pdf(paste0(Sys.Date(), "_Killfish_tissue_atlas_elbowplot.pdf"), height = 5, width= 6)
ElbowPlot(killi.atlas.clean, ndims = 50)
dev.off()

################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# While this gives us a good rough idea of the number of PCs needed to be include
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - the principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- killi.atlas.clean[["pca"]]@stdev / sum(killi.atlas.clean[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 40

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 25

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 25

# Based on these metrics, first 20 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_Killfish_tissue_atlas_elbowplot_threshold_analysis.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
################################################################################



# Run UMAP using PCA dimensionality reduction
killi.atlas.clean <- RunUMAP(killi.atlas.clean, reduction = 'pca', dims = 1:pcs)


################    summary UMAPs  ################
# By biological group
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsSinglets_UMAP_color_by_Condition.pdf", sep = "_"), height = 5, width = 5.8)
DimPlot(killi.atlas.clean, reduction = "umap", group.by = "Group", shuffle = T)
dev.off()

# By sex
my.umap.colors.sex <- c(alpha("deeppink"    , alpha = 0.3 ) ,
                        alpha("deepskyblue" , alpha = 0.3 ) )
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsSinglets_UMAP_color_by_Sex.pdf", sep = "_"), height = 5, width = 5.8)
DimPlot(killi.atlas.clean, reduction = "umap", group.by = "Sex", shuffle = T, cols = my.umap.colors.sex)
dev.off()

# By Batch
my.umap.colors.group <- c(alpha("firebrick3"    , alpha = 0.2 ) ,
                          alpha("darkturquoise" , alpha = 0.2 ) ,
                          alpha("goldenrod3"    , alpha = 0.2 ) )
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsSinglets_UMAP_color_by_Batch.pdf", sep = "_"), height = 5, width = 5.8)
DimPlot(killi.atlas.clean, reduction = "umap", group.by = "Batch", shuffle = T, cols = my.umap.colors.group)
dev.off()
### Huge batch effect!!!!
#####################################################


###############################################################################################################################################################
##### 6. Seurat RPCA integration to remove batch effects
# https://satijalab.org/seurat/articles/integration_rpca.html

# split the dataset into a list of seurat objects based on batch
killi.list <- SplitObject(killi.atlas.clean, split.by = "Batch")

# normalize and identify variable features for each dataset independently
killi.list <- lapply(X = killi.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = killi.list)

killi.list <- lapply(X = killi.list,
                     FUN = function(x) {
                       x <- ScaleData(x, features = features, verbose = FALSE)
                       x <- RunPCA(x, features = features, verbose = FALSE)
                     })

killi.anchors  <- FindIntegrationAnchors(object.list = killi.list, anchor.features = features, reduction = "rpca", k.anchor = 20)
killi.combined <- IntegrateData(anchorset = killi.anchors)

killi.combined <- ScaleData(killi.combined, verbose = FALSE)
killi.combined <- RunPCA(killi.combined, npcs = 30, verbose = FALSE)


################################################################################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# While this gives us a good rough idea of the number of PCs needed to be include
# We can calculate where the principal components start to elbow by taking the larger value of:
#    - The point where the principal components only contribute 5% of standard deviation
#    - the principal components cumulatively contribute 90% of the standard deviation.
#    - The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# Determine percent of variation associated with each PC
pct <- killi.combined[["pca"]]@stdev / sum(killi.combined[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1 # 24

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2 # 24

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs # 24

# Based on these metrics, first 20 PCs to generate the clusters.
# We can plot the elbow plot again and overlay the information determined using our metrics:

# Create a dataframe with values
plot_df <- data.frame(pct  = pct,
                      cumu = cumu,
                      rank = 1:length(pct))

# Elbow plot to visualize
pdf(paste0(Sys.Date(), "_Killfish_tissue_atlas_elbowplot_threshold_analysis_postRPCA.pdf"), height = 5, width= 6)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
################################################################################



killi.combined <- RunUMAP(killi.combined, reduction = "pca", dims = 1:pcs)

# Visualization
DimPlot(killi.combined, reduction = "umap", group.by = "Batch", shuffle = T, cols = my.umap.colors.group)


################    summary UMAPs post-batch correction ################
# By biological group
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsSinglets_UMAP_color_by_Condition_postRPCA.pdf", sep = "_"), height = 5, width = 5.8)
DimPlot(killi.combined, reduction = "umap", group.by = "Group", shuffle = T)
dev.off()

# By sex
my.umap.colors.sex <- c(alpha("deeppink"    , alpha = 0.3 ) ,
                        alpha("deepskyblue" , alpha = 0.3 ) )
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsSinglets_UMAP_color_by_Sex_postRPCA.pdf", sep = "_"), height = 5, width = 5.8)
DimPlot(killi.combined, reduction = "umap", group.by = "Sex", shuffle = T, cols = my.umap.colors.sex)
dev.off()

# By Batch
my.umap.colors.group <- c(alpha("firebrick3"    , alpha = 0.2 ) ,
                          alpha("darkturquoise" , alpha = 0.2 ) ,
                          alpha("goldenrod3"    , alpha = 0.2 ) )
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsSinglets_UMAP_color_by_Batch_postRPCA.pdf", sep = "_"), height = 5, width = 5.8)
DimPlot(killi.combined, reduction = "umap", group.by = "Batch", shuffle = T, cols = my.umap.colors.group)
dev.off()
#####################################################################


# DefaultAssay(killi.combined) <- "RNA"
# 
# FeaturePlot(object = killi.combined, features = c("LOC107372911")  , cols = c("cornsilk", "red"), order = T  )# Mpx, neutrophil/granulocyte marker
FeaturePlot(object = killi.combined, features = c("LOC107383970", # LOC107383970 Transferrin;
                                                  "fgg"         , # fibrinogen
                                                  "fgb"           # fibrinogen
) , cols = c("cornsilk", "red"), order = T )
# 
# FeaturePlot(object = killi.combined, features = c("marco")  , cols = c("cornsilk", "red"), order = T  )
# FeaturePlot(object = killi.combined, features = c("cd3e")   , cols = c("cornsilk", "red"), order = T  )
# FeaturePlot(object = killi.combined, features = c("cd79a")  , cols = c("cornsilk", "red"), order = T  )
# FeaturePlot(object = killi.combined, features = c("epor")   , cols = c("cornsilk", "red"), order = T  )

########################################################
# Seurat clustering

DefaultAssay(killi.combined) <- "integrated"

# calculate neighbors
killi.combined <- FindNeighbors(killi.combined, dims = 1:pcs)

killi.combined <- FindClusters(object = killi.combined,
                               resolution = 1,
                               n.start = 100)
# Number of communities: 51

killi.combined <- FindClusters(object = killi.combined,
                               resolution = 0.8,
                               n.start = 100)
# Number of communities: 46

killi.combined <- FindClusters(object = killi.combined,
                               resolution = 0.6,
                               n.start = 100)
# Number of communities: 42

killi.combined <- FindClusters(object = killi.combined,
                               resolution = 0.4,
                               n.start = 100)
# Number of communities: 38

killi.combined <- FindClusters(object = killi.combined,
                               resolution = 0.3,
                               n.start = 100)
# Number of communities: 31

killi.combined <- FindClusters(object = killi.combined,
                               resolution = 0.2,
                               n.start = 100)
# Number of communities: 26

killi.combined <- FindClusters(object = killi.combined,
                               resolution = 0,
                               n.start = 100)
# Number of communities: 1


# clustree: Visualise Clusterings at Different Resolutions
# https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html

pdf(paste(Sys.Date(),"Clustering_Tree_res_Singlets_CLEAN.pdf", sep = "_"), height = 15, width = 30)
clustree(killi.combined, prefix = "integrated_snn_res.")
dev.off()

save(killi.combined, file = paste(Sys.Date(),"Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_CLEAN.RData",sep = "_"))


pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsUMAP_clusters_1.0_res_Singlets_CLEAN.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.combined, label = TRUE, group.by = "integrated_snn_res.1")
dev.off()

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsUMAP_clusters_0.8_res_Singlets_CLEAN.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.combined, label = TRUE, group.by = "integrated_snn_res.0.8")
dev.off()

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsUMAP_clusters_0.6_res_Singlets_CLEAN.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.combined, label = TRUE, group.by = "integrated_snn_res.0.6")
dev.off()

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsUMAP_clusters_0.4_res_Singlets_CLEAN.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.combined, label = TRUE, group.by = "integrated_snn_res.0.4")
dev.off()

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsUMAP_clusters_0.3_res_Singlets_CLEAN.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.combined, label = TRUE, group.by = "integrated_snn_res.0.3")
dev.off()

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsUMAP_clusters_0.2_res_Singlets_CLEAN.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(killi.combined, label = TRUE, group.by = "integrated_snn_res.0.2")
dev.off()

write.table(killi.combined@meta.data, file = paste(Sys.Date(),"Seurat_PCA_SNN_clusters_and_Metadata_Singlets_CLEAN_integration.txt", sep = "_"), sep = "\t", quote = F, col.names = F, row.names = T)



################################################################################################################################################################
##### 7. Calculate cluster markers and plot potential marker gene expression

# clean and load
load('2022-10-12_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_CLEAN.RData')

# Need to re-run SCT prior to marker calculation on the integrated data [https://github.com/satijalab/seurat/issues/5321]
# After integration, the slot is compromised/missing information
killi.combined[["SCT"]] <- NULL
killi.combined     <- SCTransform(object = killi.combined, vars.to.regress = c("nFeature_RNA", "percent.mito", "Batch"), variable.features.n = 5000)
killi.combined
# An object of class Seurat 
# 38740 features across 81432 samples within 3 assays 
# Active assay: SCT (18370 features, 5000 variable features)
# 2 other assays present: RNA, integrated
# 2 dimensional reductions calculated: pca, umap


############### 0.2
# find markers for every cluster compared to all remaining cells, report only the positive ones
killi.combined <- SetIdent(object = killi.combined, value = 'integrated_snn_res.0.2')

killi.markers_0.2 <- FindAllMarkers(killi.combined, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)

# get top 15
killi.markers_0.2 %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) -> top15_0.2

write.table(data.frame(top15_0.2), file = paste(Sys.Date(),"Seurat_top15_markers_integrated_snn_res.0.2_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(killi.markers_0.2[killi.markers_0.2$p_val_adj < 0.05,], file = paste(Sys.Date(),"Seurat_ALL_FDR5_markers_integrated_snn_res.0.2_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
save(killi.markers_0.2, file = paste(Sys.Date(),"Seurat_markers_integrated_snn_res.0.2.RData", sep = "_"))

my.killi.heat_0.2 <- DoHeatmap(killi.combined, features = top15_0.2$gene, group.by = 'integrated_snn_res.0.2',   size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))  + theme(axis.text.y = element_text(size = 1))

png(paste0(Sys.Date(),"Top15_marker_heatmap_integrated_snn_res_0.2.png"), width = 40, height = 20, units = "cm", res = 300)
my.killi.heat_0.2
dev.off()

# get top 5
killi.markers_0.2 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> top5_0.2

my.killi.heat.2_0.2 <- DoHeatmap(killi.combined, features = top5_0.2$gene, group.by = 'integrated_snn_res.0.2',   size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))  + theme(axis.text.y = element_text(size = 1))

png(paste0(Sys.Date(),"Top5_marker_heatmap_integrated_integrated_res_0.2.png"), width = 40, height = 20, units = "cm", res = 300)
my.killi.heat.2_0.2
dev.off()


############### 0.3
# find markers for every cluster compared to all remaining cells, report only the positive ones
killi.combined <- SetIdent(object = killi.combined, value = 'integrated_snn_res.0.3')

killi.markers_0.3 <- FindAllMarkers(killi.combined, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)

# get top 15
killi.markers_0.3 %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) -> top15_0.3

write.table(data.frame(top15_0.3), file = paste(Sys.Date(),"Seurat_top15_markers_integrated_snn_res.0.3_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(killi.markers_0.3[killi.markers_0.3$p_val_adj < 0.05,], file = paste(Sys.Date(),"Seurat_ALL_FDR5_markers_integrated_snn_res.0.3_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
save(killi.markers_0.3, file = paste(Sys.Date(),"Seurat_markers_integrated_snn_res.0.3.RData", sep = "_"))

my.killi.heat_0.3 <- DoHeatmap(killi.combined, features = top15_0.3$gene, group.by = 'integrated_snn_res.0.3',   size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))  + theme(axis.text.y = element_text(size = 1))

png(paste0(Sys.Date(),"Top15_marker_heatmap_integrated_snn_res_0.3.png"), width = 40, height = 20, units = "cm", res = 300)
my.killi.heat_0.3
dev.off()

# get top 5
killi.markers_0.3 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> top5_0.3

my.killi.heat.2_0.3 <- DoHeatmap(killi.combined, features = top5_0.3$gene, group.by = 'integrated_snn_res.0.3',   size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))  + theme(axis.text.y = element_text(size = 1))

png(paste0(Sys.Date(),"Top5_marker_heatmap_integrated_integrated_res_0.3.png"), width = 40, height = 20, units = "cm", res = 300)
my.killi.heat.2_0.3
dev.off()


############### 0.4
# find markers for every cluster compared to all remaining cells, report only the positive ones
killi.combined <- SetIdent(object = killi.combined, value = 'integrated_snn_res.0.4')

killi.markers_0.4 <- FindAllMarkers(killi.combined, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)

# get top 15
killi.markers_0.4 %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) -> top15_0.4

write.table(data.frame(top15_0.4), file = paste(Sys.Date(),"Seurat_top15_markers_integrated_snn_res.0.4_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(killi.markers_0.4[killi.markers_0.4$p_val_adj < 0.05,], file = paste(Sys.Date(),"Seurat_ALL_FDR5_markers_integrated_snn_res.0.4_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
save(killi.markers_0.4, file = paste(Sys.Date(),"Seurat_markers_integrated_snn_res.0.4.RData", sep = "_"))

my.killi.heat_0.4 <- DoHeatmap(killi.combined, features = top15_0.4$gene, group.by = 'integrated_snn_res.0.4',   size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))  + theme(axis.text.y = element_text(size = 1))

png(paste0(Sys.Date(),"Top15_marker_heatmap_integrated_integrated_res_0.4.png"), width = 40, height = 20, units = "cm", res = 300)
my.killi.heat_0.4
dev.off()

# get top 5
killi.markers_0.4 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> top5_0.4

my.killi.heat.2_0.4 <- DoHeatmap(killi.combined, features = top5_0.4$gene, group.by = 'integrated_snn_res.0.4',   size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))  + theme(axis.text.y = element_text(size = 1))

png(paste0(Sys.Date(),"Top5_marker_heatmap_integrated_integrated_res_0.4.png"), width = 40, height = 20, units = "cm", res = 300)
my.killi.heat.2_0.4
dev.off()


############### 0.8
# find markers for every cluster compared to all remaining cells, report only the positive ones
killi.combined <- SetIdent(object = killi.combined, value = 'integrated_snn_res.0.8')

killi.markers_0.8 <- FindAllMarkers(killi.combined, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)

# get top 15
killi.markers_0.8 %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) -> top15_0.8

write.table(data.frame(top15_0.8), file = paste(Sys.Date(),"Seurat_top15_markers_integrated_snn_res.0.8_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(killi.markers_0.8[killi.markers_0.8$p_val_adj < 0.05,], file = paste(Sys.Date(),"Seurat_ALL_FDR5_markers_integrated_snn_res.0.8_Clustering.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = F)
save(killi.markers_0.8, file = paste(Sys.Date(),"Seurat_markers_integrated_snn_res.0.8.RData", sep = "_"))

my.killi.heat_0.8 <- DoHeatmap(killi.combined, features = top15_0.8$gene, group.by = 'integrated_snn_res.0.8',   size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))  + theme(axis.text.y = element_text(size = 1))

png(paste0(Sys.Date(),"Top15_marker_heatmap_integrated_integrated_res_0.8.png"), width = 40, height = 20, units = "cm", res = 300)
my.killi.heat_0.8
dev.off()

# get top 5
killi.markers_0.8 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> top5_0.8

my.killi.heat.2_0.8 <- DoHeatmap(killi.combined, features = top5_0.8$gene, group.by = 'integrated_snn_res.0.8',   size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))  + theme(axis.text.y = element_text(size = 1))

png(paste0(Sys.Date(),"Top5_marker_heatmap_integrated_integrated_res_0.8.png"), width = 40, height = 20, units = "cm", res = 300)
my.killi.heat.2_0.8
dev.off()

##########################################################################################
# requireNamespace("ape", quietly = TRUE)

killi.combined <- BuildClusterTree(object = killi.combined,
                                   dims = 1:24,
                                   assay = "integrated")


my.tree.0.8 <- Tool(object = killi.combined, slot = 'BuildClusterTree')

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsUMAP_clusters_0.8_ClusterTree.pdf", sep = "_"), height = 5, width = 10)
plot(my.tree.0.8, direction = "downwards")
dev.off()


##########################################################################################
# Resolve major cell types/clusters identity using known markers

# DefaultAssay(killi.combined) <- "RNA" # plot using the untransformed expression

##### Macrophages
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsFeature_plot_Singlets_Macrophage_markers.pdf", sep = "_"), height = 15, width = 13)
FeaturePlot(object = killi.combined, features = c("csf1r"        ,
                                                  "LOC107375663" ,  # LOC107375663 Csf1r homolog
                                                  "marco"        ,
                                                  "LOC107394065" ,  # LOC107394065: mpeg1
                                                  "LOC107387243" ,  # LOC107387243 mfap4
                                                  "lyz"             # lyzosyme
), cols = c("cornsilk", "red"), order = T )
dev.off()


##### Neutrophils
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsFeature_plot_Singlets_Neutrophil_Granulocyte_markers.pdf", sep = "_"), height = 10, width = 13)
FeaturePlot(object = killi.combined, features = c("LOC107372911", # Mpx, neutrophil/granulocyte marker
                                                  "nccrp1"      , # Neutrophil enriched https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3760823/; https://www.frontiersin.org/articles/10.3389/fimmu.2020.559555/full#h3 # Tang, JEM 2017
                                                  "LOC107386542"  # mmp13a, collagenase 3-like
                                                  
)  , cols = c("cornsilk", "red"), order = T )
dev.off()


##### B-cells
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsFeature_plot_Singlets_Bcell_markers.pdf", sep = "_"), height = 10, width = 13)
FeaturePlot(object = killi.combined, features = c("cd79a"       , #
                                                  "LOC107387651", # LOC107387651	Cd79b; cd79b CD79b molecule, immunoglobulin-associated beta
                                                  "LOC107394448", # 	immunoglobulin lambda-like polypeptide 1/igic1s1
                                                  "LOC107374380"  # 	uncharacterized protein LOC107374380/ighz
                                                  
) , cols = c("cornsilk", "red"), order = T   )
dev.off()


##### T-cells
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsFeature_plot_Singlets_Tcell_markers.pdf", sep = "_"), height = 15, width = 13)
FeaturePlot(object = killi.combined, features = c("cd3e"        ,
                                                  "LOC107372866", # LOC107372866 T-cell surface glycoprotein CD3 delta chain-like [Cd3g by blast]
                                                  "LOC107379481", # LOC107379481 cd4-1 CD4-1 molecule
                                                  "LOC107385466", # LOC107385466 cd8a CD8a molecule
                                                  "LOC107392881"  # LOC107392881 tyrosine-protein kinase ZAP-70 # Tang, JEM 2017
) , cols = c("cornsilk", "red"), order = T   )
dev.off()

##### Mast cells
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsFeature_plot_Singlets_MastCells_markers.pdf", sep = "_"), height = 10, width = 13)
FeaturePlot(object = killi.combined, features = c("LOC107377762", # LOC107377762	dicentracin-like
                                                  "LOC107377763", # LOC107377763  dicentracin-like
                                                  "LOC107377764"  # LOC107377764	pleurocidin-like peptide WF3
) , cols = c("cornsilk", "red"), order = T   )
dev.off()
# 'piscidins', a previously undiscovered family of peptide antibiotics isolated from fish, reside in mast cells,
# an immune cell of uncertain function that is present in all vertebrate classes


##### Hepatocytes
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsFeature_plot_Singlets_Hepatocytes.pdf", sep = "_"), height = 20, width = 13)
FeaturePlot(object = killi.combined, features = c("LOC107383970", # LOC107383970 Transferrin;
                                                  "ttr"         , # ttr: transthyretin
                                                  "LOC107378176", # LOC107378176: alpha-1-antitrypsin
                                                  "LOC107396791", # "LOC107396791" alpha 2-HS glycoprotein/AHSG
                                                  "fgg"         , # fibrinogen
                                                  "fgb"         , # fibrinogen
                                                  "LOC107395221", # apoa1b
                                                  "LOC107379396"  # apoa2
) , cols = c("cornsilk", "red"), order = T )
dev.off()

##### Erythrocytes
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsFeature_plot_Singlets_Erythrocytes_markers.pdf", sep = "_"), height = 10, width = 13)
FeaturePlot(object = killi.combined, features = c("hemgn"       , # hemogen regulates erythrocyte differentiation. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124579/
                                                  "epor"        , # epo receptor
                                                  "LOC107384348", # cahz carbonic anhydrase, erythrocyte https://www.nature.com/articles/s41467-017-02305-6
                                                  "LOC107390720"  # alas2 aminolevulinate, delta-, synthase 2 [erythroid-specific] https://www.nature.com/articles/s41467-017-02305-6
) , cols = c("cornsilk", "red") , order = T  )
dev.off()
################################################################################################################################################################

################################################################################################################################################################
##### 8. Try using scSorter "gates" to annotate cell types
# https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html

############ Section 1 - Generate annotation object for scSorter ############
# Prepare the annotation file for scSorter

# Read in marker file
drer.markers <- read.csv("../2022-10-10_Zebrafish_cell_markers_forAnnot.txt", sep = "\t", header = T)

# Read in homology file for killifish/zebrafish (best zebrafish hit to killifish to get conversion)
zebra.homol <- read.csv("../Zebrafish_alignment/2022-03-15_Zebrafish_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-5_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)

# create marker list object to summarize all marker genes
my.marker.list <- vector(mode = "list",length = ncol(drer.markers))
names(my.marker.list) <- colnames(drer.markers)

for (i in 1:ncol(drer.markers)) {
  
  # get unique markers from zebrafish
  my.markers <- unique(drer.markers[,i])
  my.markers <- setdiff(my.markers,"") # remove empty string
  
  # get closest Nfur corresponding gene (unique as well)
  my.marker.list[[i]] <- unique(zebra.homol$Nfur_Symbol[zebra.homol$DanRer_Symbol %in% my.markers])
  
}

names(my.marker.list)
# [1] "HSPCs"                "Thrombocytes"         "Erythrocytes"          "Neutrophils"          "Macrophages"          "B_cells"
# [7] "T_cells"              "NK_cells"         "Kidney_progenitors"   "Kidney_prox_tubule"   "Kidney_distal_tubule" "Epithelial"
# [13] "Kidney_Multiciliated" "Endothelial"         "Fibroblasts"          "Hepatocytes"          "Cholangiocytes"

# Add in manually canonical/known markers with killi genes (may be redundant)
my.marker.list$Macrophages <- unique(c(my.marker.list$Macrophages, c("csf1r"        ,
                                                                     "LOC107375663" ,  # LOC107375663 Csf1r homolog
                                                                     "marco"        ,
                                                                     "LOC107394065" ,  # LOC107394065: mpeg1
                                                                     "LOC107387243" ,  # LOC107387243 mfap4
                                                                     "lyz"             # lyzozyme
)))
my.marker.list$Neutrophils <- unique(c(my.marker.list$Neutrophils, c("LOC107372911", # Mpx, neutrophil/granulocyte marker
                                                                     "nccrp1"        # Neutrophil enriched https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3760823/; https://www.frontiersin.org/articles/10.3389/fimmu.2020.559555/full#h3 # Tang, JEM 2017
) ))
my.marker.list$Hepatocytes <- unique(c(my.marker.list$Hepatocytes, c("LOC107383970", # LOC107383970 Transferrin;
                                                                     "ttr"         , # ttr: transthyretin
                                                                     "LOC107378176", # LOC107378176: alpha-1-antitrypsin
                                                                     "LOC107396791", # "LOC107396791" alpha 2-HS glycoprotein/AHSG
                                                                     "fgg"         , # fibrinogen
                                                                     "fgb"         , # fibrinogen
                                                                     "LOC107395221", # apoa1b
                                                                     "LOC107379396"  # apoa2
                                                                     
)))
my.marker.list$Mast_cells   <- unique(c(my.marker.list$Mast_cells    , c( "LOC107377762", # LOC107377762	dicentracin-like
                                                                          "LOC107377763", # LOC107377763  dicentracin-like
                                                                          "LOC107377764"  # LOC107377764	pleurocidin-like peptide WF3
)))

my.marker.list$Erythrocytes  <- unique(c(my.marker.list$Erythrocytes    , c( "hemgn"       , # hemogen regulates erythrocyte differentiation. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124579/
                                                                             "epor"        , # epo receptor
                                                                             "LOC107384348", # cahz carbonic anhydrase, erythrocyte https://www.nature.com/articles/s41467-017-02305-6
                                                                             "LOC107390720"  # alas2 aminolevulinate, delta-, synthase 2 [erythroid-specific] https://www.nature.com/articles/s41467-017-02305-6
)))
my.marker.list$T_cells     <- unique(c(my.marker.list$T_cells    , c("cd3e"        ,
                                                                     "LOC107372866", # LOC107372866  T-cell surface glycoprotein CD3 delta chain-like [Cd3g by blast]
                                                                     "LOC107379481", # LOC107379481 cd4-1 CD4-1 molecule
                                                                     "LOC107385466", # LOC107385466 cd8a CD8a molecule
                                                                     "LOC107392881"  # LOC107392881 tyrosine-protein kinase ZAP-70 # Tang, JEM 2017,
)))


# overlap with genes detected in our data
my.detected <- rownames(killi.combined)
my.marker.list.v2 <- vector(mode = "list",length = ncol(drer.markers))

for (i in 1:length(my.marker.list)) {
  my.marker.list.v2[[i]] <- intersect(my.marker.list[[i]],my.detected)
}
names(my.marker.list.v2) <- names(my.marker.list)

# check marker lists
lapply(my.marker.list.v2, length)

# remove any definition that didn't pick up expressed/conserved genes
my.marker.list.filt <- my.marker.list.v2[lapply(my.marker.list.v2, length) > 0]

# check which cell types are remaining
names(my.marker.list.filt)
# [1] "HSPCs"                "Thrombocytes"         "Erythrocyte"          "Neutrophils"          "Macrophages"          "B_cells"
# [7] "T_cells"              "NK_cells"             "Kidney_progenitors"   "Kidney_prox_tubule"   "Kidney_distal_tubule" "Epithelial"
# [13] "Kidney_Multiciliated" "Endothelial"          "Fibroblasts"          "Hepatocytes"          "Cholangiocytes"       "Mast_cells"

### Generate anno table
# initialize
anno           <- data.frame(cbind(rep(names(my.marker.list.filt)[1],length(my.marker.list.filt[[1]])), my.marker.list.filt[[1]], rep(2,length(my.marker.list.filt[[1]]))))
colnames(anno) <- c("Type", "Marker", "Weight")

for (i in 2:length(my.marker.list.filt)) {
  anno <- rbind(anno,
                cbind("Type" = rep(names(my.marker.list.filt)[i],length(my.marker.list.filt[[i]])),
                      "Marker" = my.marker.list.filt[[i]],
                      "Weight" = rep(2,length(my.marker.list.filt[[i]]))))
  
}


# Manually increase weight for higher confidence marker genes
# e.g. Tang, JEM 2017 markers (Figure 1b)
# + canonical markers for well defined cell types
anno[anno$Marker == "cd79a"       ,]$Weight <- 6 # B-cells
anno[anno$Marker == "LOC107387651",]$Weight <- 6 # B-cells # LOC107387651	Cd79b; cd79b CD79b molecule, immunoglobulin-associated beta

anno[anno$Marker == "csf1r"       ,]$Weight <- 6 # Macrophages
anno[anno$Marker == "LOC107375663",]$Weight <- 6 # Macrophages # LOC107375663 Csf1r homolog
anno[anno$Marker == "marco"       ,]$Weight <- 6 # Macrophages
anno[anno$Marker == "LOC107394065",]$Weight <- 6 # Macrophages  # LOC107394065: mpeg1
anno[anno$Marker == "LOC107387243",]$Weight <- 6 # Macrophages  # LOC107387243 mfap4

anno[anno$Marker == "LOC107372911",]$Weight <- 6 # Neutrophil   # Mpx, neutrophil/granulocyte marker

anno[anno$Marker == "cd3e"        ,]$Weight <- 6 # T-cells
anno[anno$Marker == "LOC107372866",]$Weight <- 6 # T-cells # LOC107372866  T-cell surface glycoprotein CD3 delta chain-like [Cd3g by blast]
anno[anno$Marker == "LOC107392881",]$Weight <- 6 # T-cells # LOC107392881  tyrosine-protein kinase ZAP-70 # Tang, JEM 2017
anno[anno$Marker == "LOC107379481",]$Weight <- 6 # T-cells # LOC107379481 cd4-1 CD4-1 molecule
anno[anno$Marker == "LOC107385466",]$Weight <- 6 # T-cells # LOC107385466 cd8a CD8a molecule

anno[anno$Marker == "LOC107377762",]$Weight <- 6 # Mast cells # LOC107377762	dicentracin-like
anno[anno$Marker == "LOC107377763",]$Weight <- 6 # Mast cells # LOC107377763  dicentracin-like
anno[anno$Marker == "LOC107377764",]$Weight <- 6 # Mast cells # LOC107377764	pleurocidin-like peptide WF3
# 'piscidins', a previously undiscovered family of peptide antibiotics isolated from fish, reside in mast cells,
# an immune cell of uncertain function that is present in all vertebrate classes

anno[anno$Marker == "LOC107383970",]$Weight <- 6 # Hepatocytes # LOC107383970 Transferrin
anno[anno$Marker == "LOC107378176",]$Weight <- 6 # Hepatocytes # LOC107378176 alpha-1-antitrypsin
anno[anno$Marker == "LOC107396791",]$Weight <- 6 # Hepatocytes # LOC107396791 alpha 2-HS glycoprotein/AHSG
anno[anno$Marker == "fgg"         ,]$Weight <- 6 # Hepatocytes # fibrinogen
anno[anno$Marker == "fgb"         ,]$Weight <- 6 # Hepatocytes # fibrinogen
anno[anno$Marker == "LOC107395221",]$Weight <- 6 # Hepatocytes # apoa1b # Wang 2022
anno[anno$Marker == "LOC107379396",]$Weight <- 6 # Hepatocytes # apoa2  # Wang 2022


# make anno objects per tissue (only cell types likely to be there)
anno.liver  <- anno[anno$Type %in% c("Thrombocytes", "Erythrocytes" , "Neutrophils"   , "Macrophages" , "B_cells"         ,
                                     "T_cells"     , "NK_cells"    , "Epithelial"    , "Endothelial" ,
                                     "Fibroblasts" , "Hepatocytes" , "Cholangiocytes", "Mast_cells"),]
anno.kidney <- anno[anno$Type %in% c("HSPCs"               , "Thrombocytes", "Erythrocytes"       , "Neutrophils"       , "Macrophages"         , "B_cells"       ,
                                     "T_cells"             , "NK_cells"    , "Kidney_progenitors", "Kidney_prox_tubule", "Kidney_distal_tubule", "Epithelial"    ,
                                     "Kidney_Multiciliated", "Endothelial" , "Fibroblasts"       ,
                                     "Mast_cells"          ),]
anno.blood  <- anno[anno$Type %in% c("HSPCs"     , "Thrombocytes", "Erythrocytes" , "Neutrophils"  , "Macrophages", "B_cells",
                                     "T_cells"   , "NK_cells"    , "Mast_cells"  ),]
anno.spleen <- anno[anno$Type %in% c("HSPCs"     , "Thrombocytes", "Erythrocytes"   , "Neutrophils" , "Macrophages"       , "B_cells"      ,
                                     "T_cells"   , "NK_cells"    , "Epithelial"    , "Endothelial"    , "Fibroblasts"  ,
                                     "Mast_cells"),]


############ Section 2 - Pre-processing the data for scSorter ############

# create Seurat subset list by tissue
killi.obj.tissue <- SplitObject(killi.combined, split.by = "Tissue")

# create list of expr matrices
expr.list        <- vector(mode = "list", length = 4)
names(expr.list) <- names(killi.obj.tissue)

anno.list <- list(anno.blood, anno.kidney, anno.liver, anno.spleen)

# Get reduced expression matrix for the run
for (i in 1: length(killi.obj.tissue)) {
  # get highly variable genes and filter out genes with non-zero expression in less than 5% of total cells.
  topgenes        <- head(VariableFeatures(killi.obj.tissue[[i]]), 2500)
  expr            <- GetAssayData(killi.obj.tissue[[i]])
  topgene_filter  <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.05
  topgenes        <- topgenes[topgene_filter]
  
  # At last, we subset the preprocessed expression data.
  # Keep all genes that can be markers, even if not for cell types in that tissue
  picked_genes    <- unique(c(anno.list[[i]]$Marker, topgenes))
  print(length(picked_genes)) # 286, 797, 808, 550
  expr.list[[i]] <- expr[rownames(expr) %in% picked_genes, ]
  
  # Now, we are ready to run scSorter.
  
}

############ Section 3 - Running scSorter ############

# run on each tissue separately
rts.blood  <- scSorter(expr.list$Blood , anno.blood)
save(rts.blood, file = paste0(Sys.Date(),"_scSorter_out_Blood.RData"))

rts.liver  <- scSorter(expr.list$Liver , anno.liver)
save(rts.liver, file = paste0(Sys.Date(),"_scSorter_out_Liver.RData"))

rts.kidney <- scSorter(expr.list$Kidney, anno.kidney)
save(rts.kidney, file = paste0(Sys.Date(),"_scSorter_out_Kidney.RData"))

rts.spleen <- scSorter(expr.list$Spleen, anno.spleen)
save(rts.spleen, file = paste0(Sys.Date(),"_scSorter_out_Spleen.RData"))


######################################################## gate back ScSorter calls to Seurat object ########################################################

# get calls from CellAssign calls
killi.combined@meta.data$ScSorter_Call <- "other" # initialize

killi.combined@meta.data[colnames(expr.list$Blood ), ]$ScSorter_Call <- rts.blood$Pred_Type
killi.combined@meta.data[colnames(expr.list$Liver ), ]$ScSorter_Call <- rts.liver$Pred_Type
killi.combined@meta.data[colnames(expr.list$Kidney), ]$ScSorter_Call <- rts.kidney$Pred_Type
killi.combined@meta.data[colnames(expr.list$Spleen), ]$ScSorter_Call <- rts.spleen$Pred_Type

my.colors.cells <- c("aquamarine"    ,
                     "bisque"        ,
                     "blue"          ,
                     "blueviolet"    ,
                     "burlywood2"    ,
                     "brown3"        ,
                     "yellow"        ,
                     "chocolate1"    ,
                     "cornflowerblue",
                     "darkgoldenrod1",
                     "deeppink"      ,
                     "firebrick4"    ,
                     "forestgreen"   ,
                     "gold"          ,
                     "lawngreen"     ,
                     "lightblue1"    ,
                     "lightcoral"    ,
                     "gray83"        ,
                     "plum1"         ,
                     "slateblue1"    ,
                     "grey"          )

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsSinglets_UMAP_color_by_ScSorter_call.pdf", sep = "_"), height = 5, width = 8)
DimPlot(killi.combined, reduction = "umap", group.by = "ScSorter_Call", cols = my.colors.cells)
dev.off()

# save annotated object
save(killi.combined, file = paste0(Sys.Date(),"_Seurat_object_with_scSorter_Raw_output.RData"))


################################################################################################################################################################
##### 9. Independently try the scType algorithm (Zebrafish markers)

# https://github.com/IanevskiAleksandr/sc-type

# Load data after environment clean up
load('2022-10-04_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_CLEAN.RData')

# Run line 531-620 from the scSorter code above to get the filtered markers
my.marker.list.filt
names(my.marker.list.filt)
# [1] "HSPCs"                "Thrombocytes"         "Erythrocytes"         "Neutrophils"          "Macrophages"          "B_cells"
# [7] "T_cells"              "NK_cells"             "Kidney_progenitors"   "Kidney_prox_tubule"   "Kidney_distal_tubule" "Epithelial"
# [13] "Kidney_Multiciliated" "Endothelial"          "Fibroblasts"          "Hepatocytes"          "Cholangiocytes"       "Mast_cells"

# split Seurat by tissue of origin
killi.obj.tissue <- SplitObject(killi.combined, split.by = "Tissue")


###################################################################################################
# Construct the cell Marker object using my.marker.list.filt (instead of reading from file)

# tissueType	cellName	geneSymbolmore1 geneSymbolmore2
# Immune system	Pro-B cells	CD27,IgD,CD24,PTPRC,PAX5,CD24,CD38,CD79A,DNTT,C10orf10,VPREB1,ARPP21,CD99,IGLL1,CD9,CD79B,TCL1A,IGLL5,HLA-DQA1,HLA-DQB1,VPREB3,IGLL5
# Immune system	Macrophages	CD68,CD163,CD14,CD11b,CD206,CD80,CD86,CD16,CD64,CCL18,CD115,CD11c,CD32,HLA-DR,MRC1,MSR1,GCA	CD56

# Initialize (regardless of tissue)
cell_markers.tmp = data.frame("cellName" = names(my.marker.list.filt),
                              "geneSymbolmore1" = rep("",length(my.marker.list.filt)),
                              "geneSymbolmore2" = rep("",length(my.marker.list.filt)))
# rownames(cell_markers.tmp) <-

# grab positive markers automatically
for (i in 1:length(my.marker.list.filt)) {
  cell_markers.tmp$geneSymbolmore1[i] <- paste0(my.marker.list.filt[[cell_markers.tmp$cellName[i]]], collapse = ",")
}

# manually edit negative markers (things that shouldn't be expressed out of their cell type)
neg.markers        <- vector(mode = "list", length = length(my.marker.list.filt))
names(neg.markers) <- names(my.marker.list.filt)

neg.markers[names(neg.markers) != "Neutrophils" ] <- lapply(neg.markers[names(neg.markers) != "Neutrophils" ], function (x){append(x,"LOC107372911")}) # mpx cannot be outside of neutrophils
neg.markers[names(neg.markers) != "Hepatocytes" ] <- lapply(neg.markers[names(neg.markers) != "Hepatocytes" ], function (x){append(x,"LOC107395221")}) # LOC107395221/Apoa1b cannot be outside of hepatocytes
neg.markers[names(neg.markers) != "Hepatocytes" ] <- lapply(neg.markers[names(neg.markers) != "Hepatocytes" ], function (x){append(x,"LOC107383970")}) # Transferrin cannot be outside of hepatocytes
neg.markers[names(neg.markers) != "Hepatocytes" ] <- lapply(neg.markers[names(neg.markers) != "Hepatocytes" ], function (x){append(x,"fgb"         )}) # gb/fibrinogen cannot be outside of hepatocytes
neg.markers[names(neg.markers) != "T_cells"     ] <- lapply(neg.markers[names(neg.markers) != "T_cells"     ], function (x){append(x,"LOC107392881")}) # LOC107392881/tyrosine-protein kinase ZAP-70  cannot be outside of T-cells # Tang, JEM 2017
neg.markers[names(neg.markers) != "T_cells"     ] <- lapply(neg.markers[names(neg.markers) != "T_cells"     ], function (x){append(x,"LOC107379481")}) # LOC107379481 cd4-1 cannot be outside of T-cells
neg.markers[names(neg.markers) != "T_cells"     ] <- lapply(neg.markers[names(neg.markers) != "T_cells"     ], function (x){append(x,"LOC107385466")}) # LOC107385466 cd8a cannot be outside of T-cells
neg.markers[names(neg.markers) != "B_cells"     ] <- lapply(neg.markers[names(neg.markers) != "B_cells"     ], function (x){append(x,"cd79a"       )}) # cd79a cannot be outside of B-cells
neg.markers[names(neg.markers) != "B_cells"     ] <- lapply(neg.markers[names(neg.markers) != "B_cells"     ], function (x){append(x,"LOC107387651")}) # LOC107387651/cd79b cannot be outside of B-cells
neg.markers[names(neg.markers) != "Macrophages" ] <- lapply(neg.markers[names(neg.markers) != "Macrophages" ], function (x){append(x,"LOC107375663")}) # LOC107375663/Csf1r homolog cannot be outside of Macrophages
neg.markers[names(neg.markers) != "Erythrocytes"] <- lapply(neg.markers[names(neg.markers) != "Erythrocytes"], function (x){append(x,"hemgn"       )}) # hemgn should not be outside of erythrocytes
neg.markers[names(neg.markers) != "Erythrocytes"] <- lapply(neg.markers[names(neg.markers) != "Erythrocytes"], function (x){append(x,"epor"        )}) # epor should not be outside of erythrocytes
neg.markers[names(neg.markers) != "Erythrocytes"] <- lapply(neg.markers[names(neg.markers) != "Erythrocytes"], function (x){append(x,"LOC107384348")}) # cahz carbonic anhydrase, erythrocyte  should not be outside of erythrocytes https://www.nature.com/articles/s41467-017-02305-6
neg.markers[names(neg.markers) != "Erythrocytes"] <- lapply(neg.markers[names(neg.markers) != "Erythrocytes"], function (x){append(x,"LOC107390720")}) # alas2 aminolevulinate, delta-, synthase 2 [erythroid-specific] should not be outside of erythrocytes https://www.nature.com/articles/s41467-017-02305-6

neg.markers[names(neg.markers) != c("Macrophages","Neutrophils")] <- lapply(neg.markers[names(neg.markers) != c("Macrophages","Neutrophils")], function (x){append(x,"lyz")}) #lyz cannot be outside of macrophages and neutrophils

# populate negative markers
for (i in 1:length(my.marker.list.filt)) {
  cell_markers.tmp$geneSymbolmore2[i] <- paste0(neg.markers[[cell_markers.tmp$cellName[i]]], collapse = ",")
}

# Now create the 4 tissue configurations
liver.cell.types  <- c("Thrombocytes", "Erythrocytes" , "Neutrophils" , "Macrophages" , "B_cells" ,"T_cells" , "NK_cells" , "Epithelial" , "Endothelial" ,"Fibroblasts" , "Hepatocytes" , "Cholangiocytes", "Mast_cells")
kidney.cell.types <- c("HSPCs" , "Thrombocytes", "Erythrocytes" , "Neutrophils" , "Macrophages" , "B_cells" ,"T_cells" , "NK_cells" , "Kidney_progenitors", "Kidney_prox_tubule", "Kidney_distal_tubule", "Epithelial" ,"Kidney_Multiciliated", "Endothelial" , "Fibroblasts" ,"Mast_cells" )
blood.cell.types  <- c("HSPCs" , "Thrombocytes", "Erythrocytes" , "Neutrophils" , "Macrophages", "B_cells","T_cells" , "NK_cells" , "Mast_cells" )
spleen.cell.types <- c("HSPCs" , "Thrombocytes", "Erythrocytes" , "Neutrophils" , "Macrophages" , "B_cells" ,"T_cells" , "NK_cells" , "Epithelial" , "Endothelial" , "Fibroblasts" , "Mast_cells")

cell_markers_tissue <- rbind(cbind("tissueType" = rep("Liver" , length(liver.cell.types  )), cell_markers.tmp[cell_markers.tmp$cellName %in% liver.cell.types  ,]),
                             cbind("tissueType" = rep("Kidney", length(kidney.cell.types )), cell_markers.tmp[cell_markers.tmp$cellName %in% kidney.cell.types ,]),
                             cbind("tissueType" = rep("Blood" , length(blood.cell.types  )), cell_markers.tmp[cell_markers.tmp$cellName %in% blood.cell.types  ,]),
                             cbind("tissueType" = rep("Spleen", length(spleen.cell.types )), cell_markers.tmp[cell_markers.tmp$cellName %in% spleen.cell.types ,]) )
###################################################################################################


######  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   #####
##### Modify code from "gene_sets_prepare" of scType to marker database in R
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021

gene_sets_prepare_killi <- function(cell_markers_tissue, cell_type){
  cell_markers = cell_markers_tissue[cell_markers_tissue$tissueType == cell_type,]
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}
######  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   #####

# Source sctype function
source('/Users/benayoun/Softwares/sc-type-1.0/R/sctype_score_.R')

# prepare gene sets
gs_list.blood  <- gene_sets_prepare_killi(cell_markers_tissue, "Blood" )
gs_list.spleen <- gene_sets_prepare_killi(cell_markers_tissue, "Spleen")
gs_list.liver  <- gene_sets_prepare_killi(cell_markers_tissue, "Liver" )
gs_list.kidney <- gene_sets_prepare_killi(cell_markers_tissue, "Kidney")

# Finally, let's assign cell types to each cluster:
# get cell-type by cell matrix
es.max.blood  = sctype_score(scRNAseqData = killi.obj.tissue$Blood [["SCT"]]@scale.data, scaled = TRUE, gs = gs_list.blood$gs_positive , gs2 = gs_list.blood$gs_negative)
es.max.spleen = sctype_score(scRNAseqData = killi.obj.tissue$Spleen[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list.spleen$gs_positive, gs2 = gs_list.spleen$gs_negative)
es.max.liver  = sctype_score(scRNAseqData = killi.obj.tissue$Liver [["SCT"]]@scale.data, scaled = TRUE, gs = gs_list.liver$gs_positive , gs2 = gs_list.liver$gs_negative)
es.max.kidney = sctype_score(scRNAseqData = killi.obj.tissue$Kidney[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list.kidney$gs_positive, gs2 = gs_list.kidney$gs_negative)

# merge the results into 1 table
es.max.merged.1         <- merge(es.max.blood    , es.max.spleen  , by="row.names", all=TRUE)
es.max.merged.2         <- merge(es.max.liver    , es.max.kidney  , by="row.names", all=TRUE)
es.max.merged           <- merge(es.max.merged.1 , es.max.merged.2, by="Row.names", all=TRUE)
rownames(es.max.merged) <- es.max.merged$Row.names
es.max.merged           <- es.max.merged[,!(colnames(es.max.merged) %in% "Row.names")]

save(es.max.merged, file = paste0(Sys.Date(),"_scType_out_Matrix.RData"))

# grab highest likelihood prediction per cell (not by cluster)
killi.combined@meta.data$scType_raw = ""

for(i in 1:ncol(es.max.merged)){
  killi.combined@meta.data$scType_raw[rownames(killi.combined@meta.data) %in% colnames(es.max.merged)[i] ] = as.character(rownames(es.max.merged)[which.max(es.max.merged[,i])])
}


pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsSinglets_UMAP_color_by_ScType_call.pdf", sep = "_"), height = 5, width = 8)
DimPlot(killi.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'scType_raw', cols = my.colors.cells)
dev.off()


# save annotated object
save(killi.combined, file = paste0(Sys.Date(),"_Seurat_object_with_scSorter_scType_Raw_output.RData"))



################################################################################################################################################################
##### 10. Use cross-pecies SingleCellNet using MWS Mouse cell atlas as the model to annotate cell types

# Load data after environment clean up
load('2022-10-13_Seurat_object_with_scSorter_scType_Raw_output.RData')

# split Seurat by tissue of origin
killi.obj.tissue <- SplitObject(killi.combined, split.by = "Tissue")


###%%%%%%%%%%%%%%### Filter and prep MWS training data  ###%%%%%%%%%%%%%%###
### Load training [https://github.com/pcahan1/singleCellNet#cs_train ; 2022-10-14]
load('stList_MWS.rda')
stList_mws

# table(stList_mws$sampTab$Tissue)
# table(stList_mws$sampTab$newAnn4)

sort(unique(stList_mws$sampTab$Tissue))
# [1] "Bladder"                        "Bone_Marrow_Mesenchyme"         "Bone-Marrow"                    "Bone-Marrow_c-kit"             
# [5] "Brain"                          "Embryonic-Mesenchyme"           "Embryonic-Stem-Cell"            "Kidney"                        
# [9] "Liver"                          "Lung"                           "MammaryGland.Involution"        "MammaryGland.Lactation"        
# [13] "MammaryGland.Pregnancy"         "MammaryGland.Virgin"            "Mesenchymal-Stem-Cell-Cultured" "Muscle"                        
# [17] "Ovary"                          "Pancreas"                       "Peripheral_Blood"               "Placenta"                      
# [21] "Prostate"                       "Small-Intestine"                "Spleen"                         "Stomach"                       
# [25] "Testis"                         "Thymus"                         "Trophoblast-Stem-Cell"          "Uterus"                        

# remake object to have clean metadata slots
mws.seurat <- CreateSeuratObject(counts = stList_mws$expDat, assay = "RNA", meta.data = stList_mws$sampTab)
mws.seurat
# An object of class Seurat 
# 8058 features across 181755 samples within 1 assay 
# Active assay: RNA (8058 features, 0 variable features)

# Subset relevant tissues
mws.seurat.sub <- subset(x = mws.seurat, subset = Tissue %in% c("Bone_Marrow_Mesenchyme", "Bone-Marrow", "Bone-Marrow_c-kit", "Kidney", "Liver", "Peripheral_Blood", "Spleen") )
mws.seurat.sub
# An object of class Seurat 
# 8058 features across 52781 samples within 1 assay 
# Active assay: RNA (8058 features, 0 variable features)

# Subset on cell types
table(mws.seurat.sub@meta.data$newAnn4)

# Identify cell types with too few examples (< 30) or unlikely to be found in our tissues
my.drop.cells <- c("hematopoietic stem progenitor cell", "mast cell", "plasma cell", # 16 cells, 4 cells and 21 cells
                   "myoblast", "osteoblasts") # not likely in kidney marrow/blood/spleen/liver

# Subset cell types to remove unusable cells
mws.seurat.sub <- subset(x = mws.seurat.sub, cells = rownames(mws.seurat.sub@meta.data)[!(mws.seurat.sub@meta.data$newAnn4 %in% my.drop.cells)] )
mws.seurat.sub
# An object of class Seurat 
# 8058 features across 52208 samples within 1 assay 
# Active assay: RNA (8058 features, 0 variable features)

# Downsample to maximum 500 cells per cell type for memory efficiency
Idents(mws.seurat.sub) <- "newAnn4"
set.seed(1234567890)
mws.small <- subset(mws.seurat.sub, downsample = 500)
mws.small
# An object of class Seurat 
# 8058 features across 12264 samples within 1 assay 
# Active assay: RNA (8058 features, 0 variable features)

rm(stList_mws, mws.seurat, mws.seurat.sub)
save(mws.small, file = paste0(Sys.Date(),"_MWS_Small_for_Training.RData"))
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###


########################## TRAINING ##############################
# Load Mouse/Killifish homology table
# (best mouse hit to killifish to get conversion)
mouse.homol <- read.csv("../Mouse_alignment/2022-10-11_Mouse_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)

# Extract gene list from query data
genes.atlas.mm <- unique(mouse.homol$Mmu_Symbol[mouse.homol$Nfur_Symbol %in% rownames(killi.combined)]) # 15413

# extract info from Tabula Muris training data
MWS.sampTab       <- mws.small@meta.data
MWS.expDat        <- as.matrix(GetAssayData(mws.small, slot = "counts"))
MWS.sampTab       <- droplevels(MWS.sampTab)
MWS.sampTab$cell  <- rownames(MWS.sampTab)

# Find genes in common to the data sets and limit analysis to these
commonGenes   <- intersect(rownames(MWS.expDat), genes.atlas.mm)
length(commonGenes) # [1] 7055
MWS.expDat     <- MWS.expDat[commonGenes,]

# Split for training and assessment, and transform training data
set.seed(123456789)
stList   <- splitCommon(sampTab = MWS.sampTab, ncells = 100, dLevel = "newAnn4")
stTrain  <- stList[[1]]
expTrain <- MWS.expDat[,rownames(stTrain)]

# Train the classifier using the TMS data
class_info <- scn_train(stTrain = stTrain, expTrain = expTrain,
                        nTopGenes = 100, nTopGenePairs = 50, nRand = 50, nTrees = 1000,
                        dLevel = "newAnn4", colName_samp = "cell")
# There are 1646 top gene pairs


# ########################## TESTING ##############################
# Assessing the classifier with heldout data Apply to held out data
stTestList <- splitCommon(sampTab=stList[[2]], ncells=100, dLevel="newAnn4") # normalize validation data so that the assessment is as fair as possible
stTest     <- stTestList[[1]]
expTest    <- MWS.expDat[commonGenes,rownames(stTest)]

# predict on held out data
classRes_val_all <- scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 2)

# Assess classifier performance on held out data
tm_heldoutassessment <- assess_comm(ct_scores  = classRes_val_all,
                                    stTrain    = unique(stTrain),
                                    stQuery    = unique(stTest),
                                    dLevelSID  = "cell",
                                    classTrain = "newAnn4",
                                    classQuery = "newAnn4",
                                    nRand = 0)

pdf(paste0(Sys.Date(),"_PR_curves_SingleCellNet_MWS_model_performance.pdf") )
plot_PRs(tm_heldoutassessment)
dev.off()

pdf(paste0(Sys.Date(),"_AUPRC_curves_SingleCellNet_MWS_model_performance.pdf") )
plot_metrics(tm_heldoutassessment)
dev.off()

# Classification result heatmap
# Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
nrand          = 2
sla            = as.vector(stTest$cluster)
names(sla)     = as.vector(stTest$cell)
slaRand        = rep("rand", nrand)
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla            = append(sla, slaRand) # include in the random cells profile created

# Attribution plot
pdf(paste0(Sys.Date(),"_Classification_result_barplot_SingleCellNet_MWS_model_performance_TESTING.pdf"), height = 10, width = 15 )
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="newAnn4", sid="cell")
dev.off()
# OK

########################## QUERY ##############################

# Extract killifish data counts
my.killi.cts <- as.data.frame(GetAssayData(object = killi.combined, slot = "counts"))
my.killi.cts$GeneName <- rownames(my.killi.cts)

# merge killi/mouse based on BLAST results
my.killi.cts.ann <- merge(my.killi.cts, unique(mouse.homol[,c("Nfur_Symbol","Mmu_Symbol")]), by.x = "GeneName", by.y = "Nfur_Symbol")

# summarize based on mouse (sum paralogs if needed)
my.killi.cts.2 <- aggregate(my.killi.cts.ann[,!(colnames(my.killi.cts.ann) %in% c("GeneName","Mmu_Symbol"))], by = list(my.killi.cts.ann$Mmu_Symbol), FUN = 'sum')
rownames(my.killi.cts.2) <- my.killi.cts.2$Group.1

# Predict
classRes_killi <- scn_predict(class_info[['cnProc']], my.killi.cts.2[,-1], nrand=50)

# Classification annotation assignment
# The annotation result can be found in a column named category in the query sample table.
stkilli      <- get_cate(classRes = classRes_killi, sampTab = killi.combined@meta.data, dLevel = "Tissue", sid = "Group", nrand = 50)

table(stkilli$category, stkilli$Tissue)

# gate back to Seurat object
sum(rownames(killi.combined@meta.data) == rownames(stkilli)) # 81432
dim(killi.combined@meta.data) # 81432    27

killi.combined@meta.data$SingleCellNet_MWS          <- stkilli$category

save(stkilli, file = paste0(Sys.Date(),"_SingleCellNET_MWS_out_Matrix.RData"))

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsSinglets_UMAP_color_by_MWS_call.pdf", sep = "_"), height = 5, width = 12)
DimPlot(killi.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'SingleCellNet_MWS')
dev.off()

# save annotated object
save(killi.combined, file = paste0(Sys.Date(),"_Seurat_object_with_scSorter_scType_SingleCellNetMWS_Raw_output.RData"))


################################################################################################################################################################
##### 11. Use CellAssign [most time/memory intensive]

# Load data after environment clean up
load('2022-10-14_Seurat_object_with_scSorter_scType_SingleCellNetMWS_Raw_output.RData')

# Run lines 683-770 to reload markers in memory
names(my.marker.list.filt)
# [1] "HSPCs"                "Thrombocytes"         "Erythrocytes"         "Neutrophils"          "Macrophages"          "B_cells"             
# [7] "T_cells"              "NK_cells"             "Kidney_progenitors"   "Kidney_prox_tubule"   "Kidney_distal_tubule" "Epithelial"          
# [13] "Kidney_Multiciliated" "Endothelial"          "Fibroblasts"          "Hepatocytes"          "Cholangiocytes"       "Mast_cells"    

# Construct tissue-relevant marker gene matrices: (1) Liver cells included, (2) Kidney cells included, (3) blood cells, (4) spleen cells
my.marker.list.liver        <- my.marker.list.filt[names(my.marker.list.filt) %in% c("Thrombocytes", "Erythrocytes" , "Neutrophils"   , "Macrophages" , "B_cells"         ,
                                                                                     "T_cells"     , "NK_cells"    , "Epithelial"    , "Endothelial" ,
                                                                                     "Fibroblasts" , "Hepatocytes" , "Cholangiocytes", "Mast_cells")]

my.marker.list.kidney       <- my.marker.list.filt[names(my.marker.list.filt) %in% c("HSPCs"               , "Thrombocytes", "Erythrocytes"       , "Neutrophils"       , "Macrophages"         , "B_cells"       ,
                                                                                     "T_cells"             , "NK_cells"    , "Kidney_progenitors", "Kidney_prox_tubule", "Kidney_distal_tubule", "Epithelial"    ,
                                                                                     "Kidney_Multiciliated", "Endothelial" , "Fibroblasts"       ,
                                                                                     "Mast_cells"          )]

my.marker.list.blood <- my.marker.list.filt[names(my.marker.list.filt) %in% c("HSPCs"     , "Thrombocytes", "Erythrocytes" , "Neutrophils"  , "Macrophages", "B_cells",
                                                                              "T_cells"   , "NK_cells"    , "Mast_cells"  )]
my.marker.list.spleen <- my.marker.list.filt[names(my.marker.list.filt) %in% c("HSPCs"     , "Thrombocytes", "Erythrocytes"   , "Neutrophils" , "Macrophages"       , "B_cells"      ,
                                                                               "T_cells"   , "NK_cells"    , "Epithelial"    , "Endothelial"    , "Fibroblasts"  ,
                                                                               "Mast_cells")]
killi.marker.mat.liver  <- marker_list_to_mat(my.marker.list.liver       )
killi.marker.mat.kidney <- marker_list_to_mat(my.marker.list.kidney      )
killi.marker.mat.blood  <- marker_list_to_mat(my.marker.list.blood       )
killi.marker.mat.spleen <- marker_list_to_mat(my.marker.list.spleen      )


################################
# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
killi.sce <- as.SingleCellExperiment(killi.combined)

# object is too big to run at once - also important to only run the cell types that can be found in the sample
# we will also split and run by tissue type
my.conditions <- unique(killi.sce$Tissue)
# [1] "Blood"  "Kidney" "Liver"  "Spleen"

my.sce.list        <- vector(mode = "list",length = length(my.conditions))
names(my.sce.list) <- my.conditions

for (i in 1:length(my.conditions)) {
  # subset based on original sample
  my.sce.list[[i]] <- killi.sce[,killi.sce$Tissue %in% my.conditions[i] ]
  
  # compute size factors using the full set of genes, before subsetting to markers for input to cellassign.
  my.sce.list[[i]] <- scran::computeSumFactors(my.sce.list[[i]])
}

save(my.sce.list, file = paste0(Sys.Date(),"_SCE_normalized_object_for_CellAssign.RData"))
################################

#######################################################
# run on each sample separately in a loop (for memory and tissue specificity reasons)
my.cellassign.fit.list        <- vector(mode = "list",length = length(my.conditions))
names(my.cellassign.fit.list) <- my.conditions

my.cellassign.fit.list$Blood <- cellassign(exprs_obj        = my.sce.list$Blood[rownames(killi.marker.mat.blood),], # has to add [[1]] since we subsetted a list (we get a list of size 1)
                                           marker_gene_info = killi.marker.mat.blood                               , # has to add [[1]] since we subsetted a list (we get a list of size 1)
                                           s                = sizeFactors(my.sce.list$Blood)                       ,
                                           learning_rate    = 1e-2                                                 ,
                                           shrinkage        = TRUE                                                 ,
                                           verbose          = TRUE                                                 ,
                                           min_delta        = 2                                                    , # The minimum log fold change a marker gene must be over-expressed by in its cell type
                                           threads          = 1                                                     )

my.cellassign.fit.list$Kidney <- cellassign(exprs_obj        = my.sce.list$Kidney[rownames(killi.marker.mat.kidney),], # has to add [[1]] since we subsetted a list (we get a list of size 1)
                                            marker_gene_info = killi.marker.mat.kidney                                , # has to add [[1]] since we subsetted a list (we get a list of size 1)
                                            s                = sizeFactors(my.sce.list$Kidney)                        ,
                                            learning_rate    = 1e-2                                                   ,
                                            shrinkage        = TRUE                                                   ,
                                            verbose          = TRUE                                                   ,
                                            min_delta        = 2                                                      , # The minimum log fold change a marker gene must be over-expressed by in its cell type
                                            threads          = 1                                                      )

my.cellassign.fit.list$Liver <- cellassign(exprs_obj        = my.sce.list$Liver[rownames(killi.marker.mat.liver),], # has to add [[1]] since we subsetted a list (we get a list of size 1)
                                           marker_gene_info = killi.marker.mat.liver                             , # has to add [[1]] since we subsetted a list (we get a list of size 1)
                                           s                = sizeFactors(my.sce.list$Liver)                     ,
                                           learning_rate    = 1e-2                                               ,
                                           shrinkage        = TRUE                                               ,
                                           verbose          = TRUE                                               ,
                                           min_delta        = 2                                                  , # The minimum log fold change a marker gene must be over-expressed by in its cell type
                                           threads          = 1                                                  )

my.cellassign.fit.list$Spleen <- cellassign(exprs_obj        = my.sce.list$Spleen[rownames(killi.marker.mat.spleen),], # has to add [[1]] since we subsetted a list (we get a list of size 1)
                                            marker_gene_info = killi.marker.mat.spleen                                , # has to add [[1]] since we subsetted a list (we get a list of size 1)
                                            s                = sizeFactors(my.sce.list$Spleen)                        ,
                                            learning_rate    = 1e-2                                                   ,
                                            shrinkage        = TRUE                                                   ,
                                            verbose          = TRUE                                                   ,
                                            min_delta        = 2                                                      , # The minimum log fold change a marker gene must be over-expressed by in its cell type
                                            threads          = 1                                                      )

save(my.cellassign.fit.list, file = paste0(Sys.Date(),"_CellAssign_Fits_per_Tissue.RData"))

# get calls from CellAssign calls
killi.combined@meta.data$CellAssign_Call <- "other" # initialize
killi.combined@meta.data[colnames(my.sce.list$Blood [rownames(killi.marker.mat.blood ),]), ]$CellAssign_Call <- my.cellassign.fit.list$Blood$cell_type
killi.combined@meta.data[colnames(my.sce.list$Kidney[rownames(killi.marker.mat.kidney),]), ]$CellAssign_Call <- my.cellassign.fit.list$Kidney$cell_type
killi.combined@meta.data[colnames(my.sce.list$Liver [rownames(killi.marker.mat.liver ),]), ]$CellAssign_Call <- my.cellassign.fit.list$Liver$cell_type
killi.combined@meta.data[colnames(my.sce.list$Spleen[rownames(killi.marker.mat.spleen),]), ]$CellAssign_Call <- my.cellassign.fit.list$Spleen$cell_type

#### Plot Raw cell Assign results
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsSinglets_UMAP_color_by_CellAssign_call.pdf", sep = "_"), height = 5, width = 8)
DimPlot(killi.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'CellAssign_Call')
dev.off()

# save annotated object
save(killi.combined, file = paste0(Sys.Date(),"_Seurat_object_with_scSorter_scType_SingleCellNetMWS_cellAssign_Raw_output.RData"))



################################################################################################################################################################
##### 12. Parse cell type predictions to assign an identity to each cell cluster

load('2022-10-20_Seurat_object_with_scSorter_scType_SingleCellNetMWS_cellAssign_Raw_output.RData')

################################################################################
# Parse results from each prediction by SNN cluster group (to apply majority and decrease noise)
# we are using the 0.8 clustering resolution to improve granularity of assignment
my.cluster.annot                 <- data.frame(matrix(0, length(unique(killi.combined@meta.data$integrated_snn_res.0.8)), 19) )
rownames(my.cluster.annot)       <- paste0("Cluster_",0:(length(unique(killi.combined@meta.data$integrated_snn_res.0.8))-1))
colnames(my.cluster.annot)       <- c("integrated_snn_res.0.8",
                                      "Top_ScSorter"        , "Top_ScSorter_Perc"        ,
                                      "Second_ScSorter"     , "Second_ScSorter_Perc"     ,
                                      "Top_ScType"          , "Top_ScType_Perc"          ,
                                      "Second_ScType"       , "Second_ScType_Perc"       ,
                                      "Top_CellAssign"      , "Top_CellAssign_Perc"      ,
                                      "Second_CellAssign"   , "Second_CellAssign_Perc"   ,
                                      "Top_SingleCellNet"   , "Top_SingleCellNet_Perc"   ,
                                      "Second_SingleCellNet", "Second_SingleCellNet_Perc",                                     
                                      "Top_tissue"          ,  "Tissue_perc"             )
my.cluster.annot$integrated_snn_res.0.8 <- 0:(length(unique(killi.combined@meta.data$integrated_snn_res.0.8))-1) # initialize

# function to parse predictions
get_top2_preds <- function (clust.preds, my.colname){
  
  # tabulate by cell type called by algorithm
  tab.res <- table(meta.snn.sub[,my.colname])
  
  # Sort to get top 2
  tab.res.sort <- sort(tab.res, decreasing = T)
  
  # parse information of top 2 most frequent predictions in the cluster
  a  <- names(tab.res.sort)[1]
  b  <- 100*round(tab.res.sort[1]/sum(tab.res.sort), digits = 4)
  
  c  <- names(tab.res.sort)[2]
  d  <- 100*round(tab.res.sort[2]/sum(tab.res.sort), digits = 4)
  
  return(c(a,b,c,d))
}

#### 
for (i in 1:nrow(my.cluster.annot)) {
  
  my.snn.clust <- my.cluster.annot$integrated_snn_res.0.8[i]
  
  # subset cluster from metadata dataframe
  meta.snn.sub <- killi.combined@meta.data[killi.combined@meta.data$integrated_snn_res.0.8 == my.snn.clust,]
  
  # parse information for ScSorter
  scSorter.parse <- get_top2_preds(meta.snn.sub, "ScSorter_Call")
  my.cluster.annot[i,]$Top_ScSorter              <- scSorter.parse[1]
  my.cluster.annot[i,]$Top_ScSorter_Perc         <- scSorter.parse[2]
  my.cluster.annot[i,]$Second_ScSorter           <- scSorter.parse[3]
  my.cluster.annot[i,]$Second_ScSorter_Perc      <- scSorter.parse[4]
  
  # parse information for ScType
  scType.parse <- get_top2_preds(meta.snn.sub, "scType_raw")
  my.cluster.annot[i,]$Top_ScType                <- scType.parse[1]
  my.cluster.annot[i,]$Top_ScType_Perc           <- scType.parse[2]
  my.cluster.annot[i,]$Second_ScType             <- scType.parse[3]
  my.cluster.annot[i,]$Second_ScType_Perc        <- scType.parse[4]
  
  # parse information for CellAssign
  cellassign.parse <- get_top2_preds(meta.snn.sub, "CellAssign_Call")
  my.cluster.annot[i,]$Top_CellAssign            <- cellassign.parse[1]
  my.cluster.annot[i,]$Top_CellAssign_Perc       <- cellassign.parse[2]
  my.cluster.annot[i,]$Second_CellAssign         <- cellassign.parse[3]
  my.cluster.annot[i,]$Second_CellAssign_Perc    <- cellassign.parse[4]
  
  # parse information for SingleCellNet
  SCNET.parse <- get_top2_preds(meta.snn.sub, "SingleCellNet_MWS")
  my.cluster.annot[i,]$Top_SingleCellNet         <- SCNET.parse[1]
  my.cluster.annot[i,]$Top_SingleCellNet_Perc    <- SCNET.parse[2]
  my.cluster.annot[i,]$Second_SingleCellNet      <- SCNET.parse[3]
  my.cluster.annot[i,]$Second_SingleCellNet_Perc <- SCNET.parse[4]
  
  # tabulate by tissue
  tissue.res <- table(meta.snn.sub$Tissue)
  
  # Sort to get top
  tissue.res.sort <- sort(tissue.res, decreasing = T)
  
  # parse information
  my.cluster.annot[i,]$Top_tissue         <- names(tissue.res.sort)[1]
  my.cluster.annot[i,]$Tissue_perc        <- 100*round(tissue.res.sort[1]/sum(tissue.res.sort), digits = 4)
  
}

write.table(my.cluster.annot, file = paste(Sys.Date(),"Parsed_Cell_Annotation_Results_by_integrated_snn_res.0.8_Cluster_top.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = T)

#########################################################################################################
# populate on summarized annotation (if cluster is > 60%)
# Based on annotation of marker genes/SNN 0.8 clusters

# Automatically assign if all 3 marker based predictors agreed on top cell type
# Or if ScSorter and ScPred had same top and the first not "other" prediction for CellAssign matched
my.cluster.annot$MarkerBasedPredAgreement <-  ifelse(bitops::bitAnd(my.cluster.annot$Top_ScSorter == my.cluster.annot$Top_ScType, 
                                                                    my.cluster.annot$Top_ScSorter == ifelse(my.cluster.annot$Top_CellAssign != "other", my.cluster.annot$Top_CellAssign, my.cluster.annot$Second_CellAssign) )
                                                     >0,
                                                     my.cluster.annot$Top_ScSorter, NA)

write.table(my.cluster.annot, file = paste(Sys.Date(),"Parsed_Cell_Annotation_Results_by_integrated_snn_res.0.8_Cluster_withAutoParse.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = T)
#########################################################################################################


#########################################################################################################
# Add cell cycling info to table [to help identify progenitors]
# 2023-01-31
my.cluster.annot <- read.csv('2022-10-24_Parsed_Cell_Annotation_Results_by_integrated_snn_res.0.8_Cluster_withAutoParse.txt', sep = "\t")
my.cluster.annot$G1_percent  <- NA
my.cluster.annot$S_percent   <- NA
my.cluster.annot$G2M_percent <- NA

for (i in 1:nrow(my.cluster.annot)) {
  
  my.snn.clust <- my.cluster.annot$integrated_snn_res.0.8[i]
  
  # subset cluster from metadata dataframe
  meta.snn.sub <- killi.combined@meta.data[killi.combined@meta.data$integrated_snn_res.0.8 == my.snn.clust,]
  
  # parse information for ScSorter
  my.cluster.annot[i,]$G1_percent              <- 100*round(sum(meta.snn.sub$Phase == "G1" )/nrow(meta.snn.sub), digits = 4)
  my.cluster.annot[i,]$S_percent               <- 100*round(sum(meta.snn.sub$Phase == "S"  )/nrow(meta.snn.sub), digits = 4)
  my.cluster.annot[i,]$G2M_percent             <- 100*round(sum(meta.snn.sub$Phase == "G2M")/nrow(meta.snn.sub), digits = 4)
  
}

write.table(my.cluster.annot, file = paste(Sys.Date(),"Parsed_Cell_Annotation_Results_by_integrated_snn_res.0.8_Cluster_withAutoParse_withCellCycle.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = T)
########################################################################################################################


########################################################################################################################
# 2023-01-31/rerun for clean 2023-02-08
# Resolve ambiguous annotations
load('2022-10-20_Seurat_object_with_scSorter_scType_SingleCellNetMWS_cellAssign_Raw_output.RData')

# Read predictions and clean up object for manual curation/addition of ambiguous clusters
my.cluster.annot <- read.csv('2023-01-31_Parsed_Cell_Annotation_Results_by_integrated_snn_res.0.8_Cluster_withAutoParse_withCellCycle.txt', sep = "\t")
my.cluster.annot$MarkerBasedPredAgreement <-  ifelse(bitops::bitAnd(my.cluster.annot$Top_ScSorter == my.cluster.annot$Top_ScType, 
                                                                    my.cluster.annot$Top_ScSorter == ifelse(my.cluster.annot$Top_CellAssign != "other", my.cluster.annot$Top_CellAssign, my.cluster.annot$Second_CellAssign) ) >0,
                                                     my.cluster.annot$Top_ScSorter, NA)

# clean up table for manual additions
my.cluster.annot.v2 <- my.cluster.annot[,c("integrated_snn_res.0.8", "MarkerBasedPredAgreement")]
my.cluster.annot.v2$Final_Annotation <- my.cluster.annot.v2$MarkerBasedPredAgreement
################################################################################################


## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Clusters 6, 13, 14, 20, 27

## Clusters are overwhelmingly from liver
my.cluster.annot[my.cluster.annot.v2$integrated_snn_res.0.8 %in% c(6,13,14,20,27),]
#            Top_tissue Tissue_perc MarkerBasedPredAgreement G1_percent S_percent G2M_percent
# Cluster_6     Liver       99.94                     <NA>      52.76     30.70       16.54
# Cluster_13    Liver       99.95             Erythrocytes      52.17     33.24       14.58
# Cluster_14    Liver      100.00                     <NA>      57.58     30.97       11.45
# Cluster_20    Liver       99.75                     <NA>      77.68     13.71        8.61
# Cluster_27    Liver       98.35                     <NA>      78.56     11.84        9.60

# Prediction vary between erythrocytes and hepatocytes
# They express both erythrocytes and hepatocyte markers
clust_6.13.14.20.27 <- subset(killi.combined, subset = integrated_snn_res.0.8 %in% c("6", "13", "14", "20", "27"))
# "LOC107383970", # LOC107383970 Transferrin;
# "LOC107378176", # LOC107378176: alpha-1-antitrypsin
# "LOC107395221"  # LOC107395221/Apoa1b 
# "LOC107384348", # cahz carbonic anhydrase, erythrocyte https://www.nature.com/articles/s41467-017-02305-6
# "LOC107390720"  # alas2 aminolevulinate, delta-, synthase 2 [erythroid-specific] https://www.nature.com/articles/s41467-017-02305-6

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohorts_Cluster6_13_14_20_27_Hepa_erythro_markers.pdf", sep = "_"), height = 4, width = 6)
DotPlot(clust_6.13.14.20.27, features = c("LOC107383970", "fgg", "LOC107378176","LOC107395221",
                                          "hemgn", "epor","LOC107384348","LOC107390720"), dot.scale = 8, group.by = "integrated_snn_res.0.8") +RotatedAxis()
dev.off()

# Based on: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5790054/
# putative efferocytosing hepatocytes
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 6, ]$Final_Annotation <- "Hepatocytes_Efferocytosing"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 13,]$Final_Annotation <- "Hepatocytes_Efferocytosing"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 14,]$Final_Annotation <- "Hepatocytes_Efferocytosing"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 20,]$Final_Annotation <- "Hepatocytes_Efferocytosing"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 27,]$Final_Annotation <- "Hepatocytes_Efferocytosing"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Clusters 22
# Macrophages top choice for scType and CellAssign (82 and 94%)
clust_17.22 <- subset(killi.combined, subset = integrated_snn_res.0.8 %in% c("17", "22"))

# Weak expression of marco and mpeg1
# "LOC107394065",  # LOC107394065: mpeg1
# "LOC107387243",  # LOC107387243 mfap4
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohorts_Cluster17_22_macrophages_markers.pdf", sep = "_"), height = 4, width = 6)
DotPlot(clust_17.22, features = c("csf1r","marco", "LOC107394065", "LOC107387243"), dot.scale = 8, group.by = "integrated_snn_res.0.8") +RotatedAxis()
dev.off()

# Many complement related genes, may be more activated
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 22, ]$Final_Annotation <- "Macrophages"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Cluster 24, 49
# 24: Kidney proximal tubule top choice for scType and CellAssign (78 and 94%)
# 99% from the kidney
# 49: Kidney proximal tubule top choice for scSorter and CellAssign (97 and 99%)
# 100% from the kidney
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 24, ]$Final_Annotation <- "Kidney_prox_tubule"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 49, ]$Final_Annotation <- "Kidney_prox_tubule"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Cluster 30, 43, 46
# 30: Top predictions are Neutrophil (scSorter), Hepatocytes (ScType), Neutrophils (CellAssign) and erythroblast (singleCellNet)
# Does not express Mpx, the key neutrophil marker
# come >94% from liver tissue

# 43: Epithelial top choice for scSorter, scType and SingleCellNet
# 100% from the liver
# https://pubmed.ncbi.nlm.nih.gov/26785307/
# "The liver contains two types of epithelial cells, namely, hepatocytes and cholangiocytes"

# 46: Chologiancyte top choice for scSorter, and CellAssign; Epithelial top choice for scType and SingleCellNet
# 98% from the liver

my.cluster.annot[my.cluster.annot.v2$integrated_snn_res.0.8 %in% c(30, 43, 46),]
#                integrated_snn_res.0.8   Top_ScSorter Top_ScSorter_Perc Second_ScSorter Second_ScSorter_Perc  Top_ScType Top_ScType_Perc Second_ScType Second_ScType_Perc
# Cluster_30                     30    Neutrophils             63.33    Erythrocytes                29.56 Hepatocytes           65.59       T_cells               7.43
# Cluster_43                     43     Epithelial             91.99    Thrombocytes                 6.27  Epithelial           60.98   Macrophages              11.50
# Cluster_46                     46 Cholangiocytes             96.51    Thrombocytes                 1.74  Epithelial           39.53  Thrombocytes              17.44
#                Top_CellAssign Top_CellAssign_Perc Second_CellAssign Second_CellAssign_Perc Top_SingleCellNet Top_SingleCellNet_Perc  Second_SingleCellNet
# Cluster_30    Neutrophils               47.66      Erythrocytes                  28.92      erythroblast                  51.70       epithelial cell
# Cluster_43          other               78.75       Endothelial                   8.71   epithelial cell                  31.71 mesenchymal stem cell
# Cluster_46 Cholangiocytes               60.47        Epithelial                  37.79   epithelial cell                  77.91         megakaryocyte
#                Second_SingleCellNet_Perc Top_tissue Tissue_perc MarkerBasedPredAgreement G1_percent S_percent G2M_percent
# Cluster_30                     10.50      Liver       94.99                     <NA>      81.42     11.47        7.11
# Cluster_43                     23.69      Liver      100.00                     <NA>      85.71     10.10        4.18
# Cluster_46                      9.30      Liver       98.26                     <NA>      93.02      2.33        4.65

## rerun lines 733-754 to get cholangiocyte markers
my.marker.list$Cholangiocytes
# [1] "anxa4"        "gramd1b"      "LOC107382657" "tm4sf4"       "gadd45b"      "LOC107397367"

# Decent expression of hepatocyte markers
clust_30.43.46 <- subset(killi.combined, subset = integrated_snn_res.0.8 %in% c("30","43", "46"))
pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohorts_Cluster30_43_46_Cholangiocytes_Hepa_markers.pdf", sep = "_"), height = 4, width = 6)
DotPlot(clust_30.43.46, features = c(my.marker.list$Cholangiocytes, 
                                     "LOC107383970", "fgg", "LOC107378176","LOC107395221"), dot.scale = 8, group.by = "integrated_snn_res.0.8") +RotatedAxis()
dev.off()

my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 30, ]$Final_Annotation <- "Hepatocytes"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 43, ]$Final_Annotation <- "Hepatocytes"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 46, ]$Final_Annotation <- "Cholangiocytes"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Cluster 41
# Kidney dsital tubule top choice for scType and CellAssign (90 and 82%)
# 99% from the kidney
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 41, ]$Final_Annotation <- "Kidney_distal_tubule"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Cluster 42
# All over the place for predictions (no agreement)

# present in solid tissues, higher proportions in more vascularized kidney and spleen
table(killi.combined@meta.data$Tissue[killi.combined@meta.data$integrated_snn_res.0.8 == 42])
# Kidney  Liver Spleen 
# 123     70    100 

# Could it be vascular-related cells?
# https://academic.oup.com/eurheartj/article/36/21/1335/2293229  
# "Smooth muscle cells can be identified by several markers of their differentiation/maturation, including α-SMA, calponin, SM22, and SM-MHC isoforms (SM1 and SM2)."
# LOC107372319: calponin 2
# acta2: α-SMA
# LOC107395633:	transgelin-like	 SM22/transgelin
# https://www.ahajournals.org/doi/10.1161/ATVBAHA.120.315107
# Myh11 encodes a myosin heavy chain protein that is specifically expressed in smooth muscle cells (SMCs)
# visualize smooth muscle markers expression as dotplots 

#Fibroblast markers
# https://biomarkerres.biomedcentral.com/articles/10.1186/s40364-020-00245-w/tables/1
# Includes a-SMA (acta2)
# # # LOC107387680   # col1a1b collagen, type I, alpha 1b [ Nothobranchius furzeri (turquoise killifish) ]
# # # LOC107391613   # LOC107391613 collagen alpha-1(I) chain-like [ Nothobranchius furzeri (turquoise killifish) ]
# platelet-derived growth factor receptor alpha	pdgfra
# platelet-derived growth factor receptor beta	pdgfrb

clust_42 <- subset(killi.combined, subset = integrated_snn_res.0.8 %in% c("42"))

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohorts_Cluster42_vascular_smooth_muscle_firboblast_markers.pdf", sep = "_"), height = 4, width = 6)
DotPlot(clust_42, features = c("acta2", "LOC107372319", "LOC107395633", "myh11",
                               "col1a2","LOC107387680", "LOC107391613", "pdgfra", "pdgfrb"), group.by = "integrated_snn_res.0.8") +RotatedAxis()
dev.off()

my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 42, ]$Final_Annotation <- "Fibroblasts"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Cluster 48
# Top hits oscillate between NK cells (ScSorter, CellAssign) and T-cells (scType)
# Second best is the opposite
# 49% from the spleen
# NK and T-cell markers

# Cluster 12 is more consensually called as "T-cells", although NK is second best call

# Cluster 12 is 70% G1, and cluster 48 is 85% S
# likely that cycling is the difference
clust_12.48 <- subset(killi.combined, subset = integrated_snn_res.0.8 %in% c("12", "48"))

markers12.48 <- FindMarkers(object = killi.combined, ident.1 = "12", ident.2 = "48", group.by = 'integrated_snn_res.0.8',  logfc.threshold = 1)

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohorts_Cluster12_48_T_NK_markers.pdf", sep = "_"), height = 4, width = 6)
DotPlot(clust_12.48, features = c("cd3e","LOC107372866", "LOC107392881",
                                  "LOC107374198", # ccl38.6, NK marker in zebrafish # Tang, JEM 2017
                                  "zbtb32"      ,# NK marker in zebrafish # Tang, JEM 2017
                                  "LOC107390210"  # si:dkey-78l4.4, NK marker in zebrafish # Tang, JEM 2017
                                  
), group.by = "integrated_snn_res.0.8", scale = F) +RotatedAxis()
# "LOC107372866", # LOC107372866  T-cell surface glycoprotein CD3 delta chain-like [Cd3g by blast]
# "LOC107392881"  # LOC107392881 tyrosine-protein kinase ZAP-70 # Tang, JEM 2017
# LOC107392029 T-cell receptor beta-2 chain C region-like [ Nothobranchius furzeri (turquoise killifish) ]
# Top 10 genes identified as cluster 48 markers
# LOC107374198 C-C motif chemokine 3-like [ Nothobranchius furzeri (turquoise killifish) ]  ///// https://www.orthodb.org/?ncbi=107374198: ortholog to ccl38.6, NK marker # Tang, JEM 2017dev.off()
dev.off()

my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 12, ]$Final_Annotation <- "NK_T_cells"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 48, ]$Final_Annotation <- "NK_T_progenitor_cells"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



############## Using proliferation information to identify leukocyte progenitors ############## 

## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Cluster 10, 11, 31
# Clearly erythrocytes, but with levels of proliferation (S/G2M), suggesting a precursor state
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 10, ]$Final_Annotation <- "Erythrocyte_Progenitors"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 11, ]$Final_Annotation <- "Erythrocyte_Progenitors"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 31, ]$Final_Annotation <- "Erythrocyte_Progenitors"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Cluster 28, 35
# Clearly neutrophils, but with levels of proliferation >85% (S/G2M), suggesting a precursor state
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 28, ]$Final_Annotation <- "Neutrophil_Progenitors"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 35, ]$Final_Annotation <- "Neutrophil_Progenitors"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Cluster 26, 29
# Clearly B-cells, but with levels of proliferation >75% (S/G2M), suggesting a precursor state
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 26, ]$Final_Annotation <- "B_Cell_Progenitors"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 29, ]$Final_Annotation <- "B_Cell_Progenitors"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############## Further clustering/analysis to identify 21, 37, 38, 39 ############## 
# the purpose of the below section  is to determine which cells clusters 21, 37, 38, 39 represent

# https://traverlab.biology.ucsd.edu/wp-content/uploads/2020/09/Stachura-and-Traver.MCB_.pdf
# Figure 4: zebrafish hemopoiesis

# retain only immune related-clusters 
# to help annotate
# includes B-cells, NK/T-cells, HSPCs, macrophages and neutrophils
# B cells = clusters 9, 18, 26, 29 (progenitors 26, 29)
# NK/T cells = clusters 12, 48    (progenitors 48)
# Macrophages = clusters 17, 22 
# Neutrophils = clusters 5, 7, 50, 28, 35 (progenitors 28, 35)
# Maybe HSPC = clusters 21
# Unknown Cells = clusters 37, 38, 39
tbd.clusters <- subset(killi.combined, subset = integrated_snn_res.0.8 %in% c("9", "18", "26", "29", "12", "48", "21", "37", "38", "39","17", "22", "5", "7", "50", "28","35"))


# Cluster 21
# No agreement between methods:
# ScSorter: HSPC, ScType:Kidney_Multiciliated, cellAssign: other, SingleCellNet: erythrocyte progenitor cell
# only 21% in G1 => proliferative
# 99% in kidney, where hemopoiesis occurs (where marrow is)

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohorts_Cluster21_and_others_HSPCs_markers.pdf", sep = "_"), height = 5, width = 6)
DotPlot(tbd.clusters, features = c("zfpm1", "mef2a", "LOC107384010", "hells"), dot.scale = 8, group.by = "integrated_snn_res.0.8") +RotatedAxis()
# "LOC107384010", # tal1 T-cell acute lymphocytic leukemia 1 [ Nothobranchius furzeri (turquoise killifish) ]
dev.off()

# Clear candidate for HSPCs
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 21, ]$Final_Annotation <- "HSPCs"

############
# transfer already established annotation to help with reading and further annotation
tbd.clusters@meta.data$annot_tmp <- 'NA'
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 9 ] <- "Bcell_9"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 18] <- "Bcell_18"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 26] <- "Bcell_Pro_26"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 29] <- "Bcell_Pro_29"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 12] <- "NK_T_cell"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 48] <- "NK_T_cell_Pro"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 21] <- "HSPCs"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 17] <- "Macrophages_17"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 22] <- "Macrophages_22"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 5 ] <- "Neutrophils_5"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 7 ] <- "Neutrophils_7"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 50] <- "Neutrophils_50"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 28] <- "Neutrophils_Pro_28"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 35] <- "Neutrophils_Pro_35"

tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 37] <- "Unknown_37"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 38] <- "Unknown_38"
tbd.clusters@meta.data$annot_tmp[tbd.clusters@meta.data$integrated_snn_res.0.8 == 39] <- "Unknown_39"

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4569048/
# "Differentiated neutrophils and macrophages express the pan-leukocyte marker L-plastin (lcp1)"
# ALl subsetted clusters express it

# Clusters 37, 38, 39
# All are majority in the kidney
# Definitive hematopoietic stem cells minimally contribute to embryonic hematopoiesis
# https://www.cell.com/cell-reports/pdf/S2211-1247%2821%2901150-5.pdf
# Lymphoid-myeloid progenitors (LMPs)

markers.37 <- FindMarkers(object = tbd.clusters, ident.1 = "37", group.by = 'integrated_snn_res.0.8', only.pos = T, logfc.threshold = 0.5, min.pct = 0.75)
head(markers.37[order(markers.37$avg_logFC, decreasing = T),], n = 10L)
#                       p_val avg_logFC pct.1 pct.2     p_val_adj
# sox4         1.267421e-192 1.0744044 0.829 0.223 2.328253e-188
# LOC107394792 2.461848e-153 1.0455148 0.945 0.403 4.522414e-149 # fabp3 fatty acid binding protein 3, muscle and heart [ Nothobranchius furzeri (turquoise killifish) ]
# csf1r         0.000000e+00 0.9574236 0.780 0.112  0.000000e+00 
# LOC107379844  0.000000e+00 0.9293865 0.815 0.097  0.000000e+00 # LOC107379844 zinc finger protein 36, C3H1 type-like 1 [ Nothobranchius furzeri (turquoise killifish) ]
# hmga2         0.000000e+00 0.9104695 0.763 0.065  0.000000e+00 
# eef1g        7.111566e-131 0.8075050 0.997 0.962 1.306395e-126 
# hsp90ab1     8.351711e-157 0.7967385 1.000 0.998 1.534209e-152 
# LOC107393303 1.032088e-114 0.7492474 0.994 0.912 1.895945e-110 # elongation factor 2-like
# LOC107384299 1.298350e-111 0.7253631 0.884 0.417 2.385069e-107 # cebpb CCAAT enhancer binding protein beta
# LOC107379384  2.648099e-85 0.6910181 0.789 0.377  4.864557e-81 # cited4b Cbp/p300-interacting transactivator, with Glu/Asp-rich carboxy-terminal domain, 4b

markers.38 <- FindMarkers(object = tbd.clusters, ident.1 = "38", group.by = 'integrated_snn_res.0.8', only.pos = T, logfc.threshold = 0.5, min.pct = 0.75)
head(markers.38[order(markers.38$avg_logFC, decreasing = T),], n = 10L)
#                      p_val avg_logFC pct.1 pct.2     p_val_adj
# LOC107386128  0.000000e+00 2.9044201 0.874 0.064  0.000000e+00 # LOC107386128 actin, cytoplasmic 1-like [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107375851 6.662420e-153 1.1661233 0.821 0.270 1.223886e-148 # cd74b CD74 molecule, major histocompatibility complex, class II invariant chain b 
# cd74         3.401093e-116 0.9595992 0.894 0.371 6.247808e-112 # 
# sftpb        1.895149e-110 0.8537407 0.924 0.377 3.481389e-106 # 
# LOC107384938  6.371884e-59 0.7841205 0.829 0.744  1.170515e-54 # LOC107384938 actin-related protein 2/3 complex subunit 4-like 
# LOC107372638  5.731196e-83 0.7620119 0.859 0.367  1.052821e-78 # H-2 class II histocompatibility antigen, E-S beta chain-like
# LOC107373955  1.774049e-45 0.7026193 0.935 0.787  3.258929e-41 # profilin-2-like
# actb          1.569557e-54 0.6095000 1.000 0.999  2.883277e-50 # 
# LOC107389947  1.033056e-35 0.5621204 0.988 0.881  1.897724e-31 # uncharacterized
# lcp1          2.354296e-29 0.5168300 0.865 0.776  4.324842e-25 # 

markers.39 <- FindMarkers(object = tbd.clusters, ident.1 = "39", group.by = 'integrated_snn_res.0.8', only.pos = T, logfc.threshold = 0.5, min.pct = 0.75)
head(markers.39[order(markers.39$avg_logFC, decreasing = T),], n = 10L)
#                       p_val avg_logFC pct.1 pct.2     p_val_adj
# LOC107379864   0.000000e+00  2.843693 0.994 0.148  0.000000e+00 # rnaset2l ribonuclease T2, like [ Nothobranchius furzeri (turquoise killifish) ]
# sftpb         1.084206e-238  2.050250 1.000 0.377 1.991687e-234 # 
# LOC107375851  7.462106e-226  1.926821 0.952 0.269 1.370789e-221 # cd74b CD74 molecule, major histocompatibility complex, class II invariant chain b
# csgr14hxorf21  0.000000e+00  1.823083 0.933 0.147  0.000000e+00 # csgr14hxorf21 chromosome sgr14 open reading frame, human CXorf21
# cd74          2.110786e-196  1.673914 0.997 0.370 3.877515e-192 # 
# LOC107379614   0.000000e+00  1.615898 0.895 0.006  0.000000e+00 # si:ch211-105c13.3 uncharacterized LOC107379614 [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107385064   0.000000e+00  1.507796 0.952 0.122  0.000000e+00 # uncharacterized
# LOC107385077  8.456807e-278  1.358694 0.984 0.220 1.553515e-273 # uncharacterized
# LOC107372638  4.965466e-167  1.330295 0.990 0.366 9.121562e-163 # H-2 class II histocompatibility antigen, E-S beta chain-like
# csf1r          0.000000e+00  1.303641 0.920 0.111  0.000000e+00 # 
##
# lyz             0  1.289568 0.755 0.086         0



pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohorts_Cluster21_37_38_39_and_related_HSPCs_and_otherblood_markers.pdf", sep = "_"), height = 6, width = 7)
DotPlot(tbd.clusters, features = c("zfpm1"       , # HSPC markers
                                   "hells"       , # HSPC markers

                                   "cd74"        , # CD74 is a regulator of hematopoietic stem cell maintenance
                                   "flt3"        , # Progenitors
                                   "sox4"        , # B-cell progenitor commitment marker
                                   "ikzf1"       , # Ikaros lymphoid

                                   "csf1r"       , # Macrophage markers
                                   "lyz"         , # Macrophage markers

                                   "LOC107372911", # Mpx, neutrophil/granulocyte marker
                                   "cd3e"        , # T-cell marker
                                   "cd79a"         # B-cell marker
                                    ),
        dot.scale = 8, group.by = "annot_tmp") +RotatedAxis()
dev.off()


### all have high cd74
# CD74 is a regulator of hematopoietic stem cell maintenance
# https://pubmed.ncbi.nlm.nih.gov/33661886/

### Sox4 (highest in 37)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3390953/
# "The Sox4 transcription factor was identified by our screen as a gene that cooperates with CREB in myeloid leukemogenesis;
# These results indicate that Sox4 and CREB cooperate and contribute to increased proliferation of hematopoietic progenitor cells."
#####
# https://ashpublications.org/blood/article/124/21/1577/88774/Sox4-Is-Required-for-the-Formation-and-Maintenance
# "study reveals a novel role for Sox4 gene in early hematopoiesis and brings important insights into the regulatory mechanisms underlying the commitment of HSCs toward multipotent progenitors."
#### 
# https://pubmed.ncbi.nlm.nih.gov/23345330/
# "Sox4 is required for the survival of pro-B cells"

### Flt3
# https://ashpublications.org/blood/article/100/5/1532/106333/The-roles-of-FLT3-in-hematopoiesis-and-leukemia
# FLT3 is expressed in a variety of human and murine cell lines of both myeloid and B-lymphoid lineage.
# In normal bone marrow, expression appears to be restricted to early progenitors
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4089881/
# Flt3 signaling regulates the proliferation, survival, and maintenance of multipotent hematopoietic progenitors that generate B cell precursors
# MPP

### csf1r
# can be expressed in multipotent progenitors
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4031967/
# Low, HSC-equivalent levels of expression of Csf1r mRNA are found on common myeloid (CMP) and common lymphoid (CLP) progenitor cells

### Ikaros/ikzf1
# https://pubmed.ncbi.nlm.nih.gov/11748838/
# Ikaros expression as a marker for lymphoid progenitors during zebrafish development

### ccr9a
# XP_015810258.1	C-C chemokine receptor type 9-like	LOC107382569	ccr9a	ENSDARP00000071978	ENSDART00000077511
# "Ikzf1 regulates embryonic T lymphopoiesis via Ccr9 and Irf4 in zebrafishIkzf1 regulates T lymphopoiesis via Ccr9 and Irf4"
# https://www.sciencedirect.com/science/article/pii/S0021925820304701

### "LOC107378674", # AIF-1
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5943248/
#   AIF loss deregulates hematopoiesis and reveals different adaptive metabolic responses in bone marrow cells and thymocytes

#### Based on all this information:
# 37 with hells (stem), Cd74 (progenitor maintenance), high sox4 (progenitor/B), high csf1r but NO lyz, Ikaros and ccr91 => likely CLP/Lymphoid progenitor
# 38 with high Cd74 (progenitor maintenance), flt3, AIF1, LOTs of actin related genes as top markers => multipotent progenitor
# 39 with high csf1r, high lyz, low Mpx (but highest of non-neutrophils) => likely myeloid progenitor

my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 37, ]$Final_Annotation <- "Lymphoid_progenitors"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 38, ]$Final_Annotation <- "Multipotent_progenitors"
my.cluster.annot.v2[my.cluster.annot.v2$integrated_snn_res.0.8 == 39, ]$Final_Annotation <- "Myeloid_progenitors"
## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

write.table(my.cluster.annot.v2, file = paste(Sys.Date(),"Parsed_Cell_Annotation_Results_by_integrated_snn_res.0.8_FINAL_CURATED.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = T)


#################################################################################################################################################################################
# Transfer annotation to Seurat object
killi.combined@meta.data$Annotation_v1 <- NA

for (i in 1:nrow(my.cluster.annot.v2)) {
  killi.combined@meta.data$Annotation_v1[killi.combined@meta.data$integrated_snn_res.0.8 ==  my.cluster.annot.v2$integrated_snn_res.0.8[i]] <- my.cluster.annot.v2$Final_Annotation[i]
}


pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsUMAP_Annotation_v1_plot.pdf", sep = "_"), height = 5, width = 9)
DimPlot(killi.combined, label = TRUE, group.by = "Annotation_v1", repel = T) + NoLegend()
dev.off()

# ==> Need to filter cells so that we remove cells that are clustering with an identity impossible based on tissue of origin
table(killi.combined@meta.data$Annotation_v1,killi.combined@meta.data$Tissue)
#                               Blood Kidney Liver Spleen
# B_Cell_Progenitors               0   1372    10     26
# B_cells                         50   2394    69   1264
# Cholangiocytes                   0      0   169      3
# Endothelial                      0    137  1348     97
# Erythrocyte_Progenitors         30   4631    15    168
# Erythrocytes                 15050   2744  3466  14191
# Fibroblasts                      0    123    70    100
# Hepatocytes                      0      1  6577     31
# Hepatocytes_Efferocytosing       4      3  8714     10
# HSPCs                            0   1164     1      6
# Kidney_distal_tubule             0    305     0      0
# Kidney_prox_tubule               0   1020     0      1
# Lymphoid_progenitors             1    331     6      8
# Macrophages                     16    865   189   1306
# Mast_cells                       2    159   122    944
# Multipotent_progenitors          0    198     5    137
# Myeloid_progenitors              0    194    16    104
# Neutrophil_Progenitors           5   1023    34     11
# Neutrophils                     31   6819   273    254
# NK_T_cells                       8   1110   148    637
# NK_T_progenitor_cells            1     53    25     76
# Thrombocytes                    81     99   123    654

killi.combined@meta.data$Tissue_QC <- "Yes"
killi.combined@meta.data$Tissue_QC[bitAnd(killi.combined@meta.data$Tissue != "Liver" , killi.combined@meta.data$Annotation_v1 == "Cholangiocytes"            )>0] <- "No"
killi.combined@meta.data$Tissue_QC[bitAnd(killi.combined@meta.data$Tissue != "Liver" , killi.combined@meta.data$Annotation_v1 == "Hepatocytes"               )>0] <- "No"
killi.combined@meta.data$Tissue_QC[bitAnd(killi.combined@meta.data$Tissue != "Liver" , killi.combined@meta.data$Annotation_v1 == "Hepatocytes_Efferocytosing")>0] <- "No"
killi.combined@meta.data$Tissue_QC[bitAnd(killi.combined@meta.data$Tissue != "Kidney", killi.combined@meta.data$Annotation_v1 == "HSPCs"                     )>0] <- "No"
killi.combined@meta.data$Tissue_QC[bitAnd(killi.combined@meta.data$Tissue != "Kidney", killi.combined@meta.data$Annotation_v1 == "Kidney_prox_tubule"        )>0] <- "No"
killi.combined@meta.data$Tissue_QC[bitAnd(killi.combined@meta.data$Tissue != "Kidney", killi.combined@meta.data$Annotation_v1 == "Lymphoid_progenitors"      )>0] <- "No"


### Keep only QC cells
killi.combined.filt <- subset(killi.combined, subset = Tissue_QC %in% "Yes")
killi.combined.filt
# An object of class Seurat 
# 38740 features across 81357 samples within 3 assays 
# Active assay: SCT (18370 features, 5000 variable features)
# 2 other assays present: RNA, integrated
# 2 dimensional reductions calculated: pca, umap

table(killi.combined.filt@meta.data$Annotation_v1,killi.combined.filt@meta.data$Tissue)
#                               Blood Kidney Liver Spleen
#     B_Cell_Progenitors             0   1372    10     26
#     B_cells                       50   2394    69   1264
#     Cholangiocytes                 0      0   169      0
#     Endothelial                    0    137  1348     97
#     Erythrocyte_Progenitors       30   4631    15    168
#     Erythrocytes               15050   2744  3466  14191
#     Fibroblasts                    0    123    70    100
#     Hepatocytes                    0      0  6577      0
#     Hepatocytes_Efferocytosing     0      0  8714      0
#     HSPCs                          0   1164     0      0
#     Kidney_distal_tubule           0    305     0      0
#     Kidney_prox_tubule             0   1020     0      0
#     Lymphoid_progenitors           0    331     0      0
#     Macrophages                   16    865   189   1306
#     Mast_cells                     2    159   122    944
#     Multipotent_progenitors        0    198     5    137
#     Myeloid_progenitors            0    194    16    104
#     Neutrophil_Progenitors         5   1023    34     11
#     Neutrophils                   31   6819   273    254
#     NK_T_cells                     8   1110   148    637
#     NK_T_progenitor_cells          1     53    25     76
#     Thrombocytes                  81     99   123    654

save(killi.combined.filt, file = paste0(Sys.Date(),"_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData"))


# get UMAP
install.packages("Polychrome")
library(Polychrome)

# create your own color palette based on `seedcolors`
P22 = createPalette(22,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P22)

pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsUMAP_Annotation_v1_plot_CLEAN.pdf", sep = "_"), height = 5, width = 6)
DimPlot(killi.combined, label = TRUE, group.by = "Annotation_v1", repel = T, cols = as.vector(P22)) + NoLegend()
dev.off()


pdf(paste(Sys.Date(),"Killifish_Tissues_3_cohortsUMAP_Annotation_v1_plot_CLEAN_withlegend.pdf", sep = "_"), height = 5, width = 10)
DimPlot(killi.combined, label = F, group.by = "Annotation_v1", repel = T, cols = as.vector(P22))
dev.off()
################################################################################################################################################################


################################################################################################################################################################
##### 13. get cell type proportions for each library

# get cell type proportions
my.freq.table <- prop.table(x = table(killi.combined.filt@meta.data$Annotation_v1, killi.combined.filt@meta.data$Condition), margin = 2)
write.table(t(my.freq.table), file = paste(Sys.Date(),"Cell_frequencies_per_sample.txt", sep = "_"), sep = "\t", quote = F, col.names = T, row.names = T)

# make beeswarm representation
my.freqs <- as.data.frame(as.matrix(t(as.matrix(my.freq.table))))
my.freqs

### Make frequency beeswarms
my.freqs$Tissue <- sapply(strsplit(as.character(my.freqs$Var1), "_"),'[',2)
my.freqs$Sex    <- sapply(strsplit(as.character(my.freqs$Var1), "_"),'[',1)

#### Blood
blood.freqs <- my.freqs[my.freqs$Tissue == "Blood",]

pdf(paste(Sys.Date(),"Killifish_Tissues_BEESWARM_Blood_Cell_Frequencies.pdf", sep = "_"), height = 5, width = 8)
beeswarm::beeswarm(Freq ~ Var2, data = blood.freqs, pc = 16, pwcol = ifelse(blood.freqs$Sex == "Male","deepskyblue","deeppink"),
                   xlab = "", ylab = "Fraction of cells in 10xGenomics single cell library", main = "Blood", las = 2, corral = "random")
dev.off()

#### Spleen
spleen.freqs <- my.freqs[my.freqs$Tissue == "Spleen",]

pdf(paste(Sys.Date(),"Killifish_Tissues_BEESWARM_Spleen_Cell_Frequencies.pdf", sep = "_"), height = 5, width = 8)
beeswarm::beeswarm(Freq ~ Var2, data = spleen.freqs, pc = 16, pwcol = ifelse(spleen.freqs$Sex == "Male","deepskyblue","deeppink"),
                   xlab = "", ylab = "Fraction of cells in 10xGenomics single cell library", main = "Spleen", las = 2, corral = "random")
dev.off()

#### Liver
Liver.freqs <- my.freqs[my.freqs$Tissue == "Liver",]

pdf(paste(Sys.Date(),"Killifish_Tissues_BEESWARM_Liver_Cell_Frequencies.pdf", sep = "_"), height = 5, width = 8)
beeswarm::beeswarm(Freq ~ Var2, data = Liver.freqs, pc = 16, pwcol = ifelse(Liver.freqs$Sex == "Male","deepskyblue","deeppink"),
                   xlab = "", ylab = "Fraction of cells in 10xGenomics single cell library", main = "Liver", las = 2, corral = "random")
dev.off()

#### Kidney
Kidney.freqs <- my.freqs[my.freqs$Tissue == "Kidney",]

pdf(paste(Sys.Date(),"Killifish_Tissues_BEESWARM_Kidney_Cell_Frequencies.pdf", sep = "_"), height = 5, width = 8)
beeswarm::beeswarm(Freq ~ Var2, data = Kidney.freqs, pc = 16, pwcol = ifelse(Kidney.freqs$Sex == "Male","deepskyblue","deeppink"),
                   xlab = "", ylab = "Fraction of cells in 10xGenomics single cell library", main = "Kidney", las = 2, corral = "random")
dev.off()


################################################################################################################################################################
##### 14. Show cell type marker expression

# Cells
table(killi.combined.filt@meta.data$Annotation_v1)
#           B_Cell_Progenitors                      B_cells               Cholangiocytes                  Endothelial 
#                         1408                         3777                          169                         1582 
#      Erythrocyte_Progenitors                 Erythrocytes                  Hepatocytes   Hepatocytes_Efferocytosing 
#                         4844                        35451                         6577                         8714 
#                        HSPCs         Kidney_distal_tubule           Kidney_prox_tubule         Lymphoid_progenitors 
#                         1164                          305                         1020                          331 
#                  Macrophages                   Mast_cells      Multipotent_progenitors          Myeloid_progenitors 
#                         2376                         1227                          340                          314 
#       Neutrophil_Progenitors                  Neutrophils                   NK_T_cells        NK_T_progenitor_cells 
#                         1073                         7377                         1903                          155 
#                 Thrombocytes Vascular_Smooth_Muscle_Cells 
#                          957                          293 

sort(unique(killi.combined.filt@meta.data$Annotation_v1))
# B_Cell_Progenitors
# B_cells
# Cholangiocytes
# Endothelial
# Erythrocyte_Progenitors
# Erythrocytes
# Fibroblasts
# Hepatocytes
# Hepatocytes_Efferocytosing
# HSPCs
# Kidney_distal_tubule
# Kidney_prox_tubule
# Lymphoid_progenitors
# Macrophages
# Mast_cells
# Multipotent_progenitors
# Myeloid_progenitors
# Neutrophil_Progenitors
# Neutrophils
# NK_T_cells
# NK_T_progenitor_cells
# Thrombocytes


# 2-3 markers per cell type
my.markers <- c("cd79a"       , #
                "LOC107387651", # LOC107387651	Cd79b; cd79b CD79b molecule
                "pax5"        , # Pax5 functions as a B cell identity factor throughout B cell development
                
                # Cholangiocyte
                "anxa4"       , # Cholangiocyte markers
                "tm4sf4"      , # Cholangiocyte markers
                "LOC107382657", # Cholangiocyte markers/klf6a Kruppel-like factor 6a

                # Endothelial
                "cpn1"        , #
                "LOC107378530", # kdrl kinase insert domain receptor like [ Nothobranchius furzeri (turquoise killifish) ]
                "LOC107375902", # cdh5 cadherin 5 [ Nothobranchius furzeri (turquoise killifish) ]
                
                # Erythrocytes
                "LOC107390720", # alas2 aminolevulinate, delta-, synthase 2 [ Nothobranchius furzeri (turquoise killifish) ]
                "LOC107384348", # cahz carbonic anhydrase, erythrocyte https://www.nature.com/articles/s41467-017-02305-6
                "LOC107378372", # hemoglobin subunit beta-A-like [ Nothobranchius furzeri (turquoise killifish) ]

                # Fibroblasts
                "acta2"       ,
                "col1a2"      ,
                "pdgfrb"      ,
                
                # Hepatocytes
                "LOC107383970" , # LOC107383970 Transferrin;
                "fgg"          , # fibrinogen gamma chain
                "LOC107378176" , # LOC107378176: alpha-1-antitrypsin

                # HSPCs
                "hells"       ,
                "zfpm1"       ,
                "LOC107384010", # tal1 T-cell acute lymphocytic leukemia 1 [ Nothobranchius furzeri (turquoise killifish) ]

                # Kidney_distal_tubule
                "LOC107375383", # cytochrome c oxidase subunit 5B, mitochondrial-like
                "LOC107377205", # sodium/potassium-transporting ATPase subunit alpha-1
                "bin2"        , # bridging integrator 2

                # Kidney_prox_tubule
                "LOC107374907", # fbp1-like fructose-1,6-bisphosphatase 1-like
                "LOC107377841", # slc22a2-like solute carrier family 22 member 2-like
                "pdzk1"       ,

                # lymphoid progenitors
                "cd74"        , # CD74 is a regulator of hematopoietic stem cell maintenance
                "sox4"        , #  B-cell progenitor commitment marker
                "ikzf1"       , # Ikaros lymphoid
                
                # Macrophages
                "csf1r"       , # Macrophages
                "marco"       , # Macrophages
                "LOC107394065", # LOC107394065: mpeg1

                # Mast Cells
                "LOC107377762", # LOC107377762	dicentracin-like
                "LOC107377763", # LOC107377763  dicentracin-like
                "LOC107377764", # LOC107377764	pleurocidin-like peptide WF3
                
                # MPP
                "flt3"        , # Progenitors
                
                # Myeloid_progenitors
                "lyz"         , # Macrophage/myeloid markers
                
                # Neutrophils
                "LOC107372911", # Mpx, neutrophil/granulocyte marker
                "LOC107386182", # il1b interleukin 1, beta [ Nothobranchius furzeri (turquoise killifish) ]
                "nccrp1"      , # Neutrophil enriched https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3760823/; https://www.frontiersin.org/articles/10.3389/fimmu.2020.559555/full#h3 # Tang, JEM 2017
                
                # T cell
                "cd3e"        ,
                "LOC107372866", # LOC107372866  T-cell surface glycoprotein CD3 delta chain-like [Cd3g by blast]
                "LOC107392881", # LOC107392881 tyrosine-protein kinase ZAP-70 # Tang, JEM 2017

                # Thrombocyte
                "LOC107377820", # thbs1b thrombospondin 1b [ Nothobranchius furzeri (turquoise killifish) ]
                "fgl2"        , # fibrinogen like 2
                "LOC107391531", # flna filamin A, alpha (actin binding protein 280) [ Nothobranchius furzeri (turquoise killifish) ]

                # Proliferation
                "pcna",
                "mki67"
                )

pdf(paste(Sys.Date(),"Killifish_Tissues_Dotplot_Cell_type_markers_CLEAN.pdf", sep = "_"), height = 8, width = 16)
DotPlot(killi.combined.filt, features = my.markers, dot.scale = 8, group.by = "Annotation_v1", cols = c("gray96","blue"))+RotatedAxis()
dev.off()


#### for heatmap, needs more balance: subset 200 cells per cell type to plot (otherwise, readability is terrible)
# https://satijalab.org/seurat/v3.0/multimodal_vignette.html
killi.combined.filt <- SetIdent(object = killi.combined.filt, value = 'Annotation_v1')
table(killi.combined.filt@meta.data$Annotation_v1)

killi.small.plot <- subset(killi.combined.filt, downsample = 200)
DefaultAssay(killi.small.plot) <- "RNA"


# Get normalized RNA values
killi.small.plot <- NormalizeData(killi.small.plot)
killi.small.plot <- FindVariableFeatures(killi.small.plot, selection.method = "vst", nfeatures = 5000)
killi.small.plot <- ScaleData(killi.small.plot, features = rownames(killi.small.plot))

pdf(paste(Sys.Date(),"Killifish_Tissues_Heatmap_Cell_type_markers_Annov2_200cells_per_Type_yellow_purple.pdf", sep = "_"), width = 17, height = 11)
DoHeatmap(killi.small.plot, features = my.markers, group.by = "Annotation_v1", group.colors  = as.vector(P22)) + theme(axis.text.y = element_text(size = 12))
dev.off()

pdf(paste(Sys.Date(),"Killifish_Tissues_Heatmap_Cell_type_markers_Annov2_200cells_per_Type.pdf", sep = "_"), width = 17, height = 11)
DoHeatmap(killi.small.plot, features = my.markers, group.by = "Annotation_v1") + theme(axis.text.y = element_text(size = 12)) + scale_fill_gradientn(colors = c("darkblue", "white", "red"))
dev.off()

# #################################################################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_Killi_Atlas_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()
