# setwd('/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2022-02-01_Killifish_scRNAseq_ATLAS_analyses/Differential_Expression/')
setwd('/Users/berenice/Dropbox/Manuscripts_and_Publications/2023/2023_Killifish_single_cell_Atlas_Bryan/TEST/edgeR_LRT')
options(stringsAsFactors = F)

### R version 4.1.2 (2021-11-01)

library('Seurat')         # SeuratObject_4.1.3
library('muscat')         # muscat_1.8.2    
library(sctransform)      # sctransform_0.3.5 
library("singleCellTK")   # singleCellTK_2.4.0
library('DESeq2')         # DESeq2_1.34.0  
library('RUVSeq')         # RUVSeq_1.28.0  
library(clusterProfiler)  # clusterProfiler_4.2.2
library(org.Dr.eg.db)     # org.Dr.eg.db_3.14.0 
require(DOSE)             # DOSE_3.20.1     
library(ggplot2)          # ggplot2_3.4.0    
library(scales)           # scales_1.2.1 

library("ComplexHeatmap") # ComplexHeatmap_2.10.0
library("bitops")         # bitops_1.0-7 
library(circlize)         # circlize_0.4.15
library(Polychrome)       # Polychrome_1.5.1

theme_set(theme_bw())   

library(edgeR)            # edgeR_3.36.0

# 2023-02-10
# Run on clean new annotation

# 2023-06-08
# Run edgeR-LRT as alternative pseudobulk method

###############################################################################################
# 0. preprocess Seurat object for use with muscat

# Load annotated cell data
load('2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData')
killi.combined.filt
# An object of class Seurat 
# 38740 features across 81357 samples within 3 assays 
# Active assay: SCT (18370 features, 5000 variable features)
# 2 other assays present: RNA, integrated
# 2 dimensional reductions calculated: pca, umap

killi.combined.filt@meta.data$Fish <- paste0(killi.combined.filt@meta.data$Sex,"_",killi.combined.filt@meta.data$Batch)

# bring RNA as main assay again
DefaultAssay(killi.combined.filt) <- "RNA"
killi.combined.filt[['SCT']] <- NULL

# convert to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
killi.atlas.sce <- as.SingleCellExperiment(killi.combined.filt)

rm(killi.combined.filt) # free some memory

###############################################################################################
# 1. Run muscat for pseudobulking and extraction of samples
# Data preparation
#### param names changed in new muscat package version (2023-06-08)
killi.atlas.sce.cl <- prepSCE(killi.atlas.sce, 
                              kid   = "Annotation_v1",  # population assignments
                              gid   = "Sex"          ,  # group IDs (ctrl/stim)
                              sid   = "Fish"         ,  # sample IDs (ctrl/stim.1234)
                              drop  = TRUE)             # drop all other colData columns


# store cluster and sample IDs, as well as the number of clusters and samples into the following simple variables:
nk  <- length(kids <- levels(killi.atlas.sce.cl$cluster_id))
ns  <- length(sids <- levels(killi.atlas.sce.cl$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(killi.atlas.sce.cl$cluster_id, killi.atlas.sce.cl$sample_id))
#             B_Cell_Progenitors B_cells Cholangiocytes Endothelial Erythrocyte_Progenitors Erythrocytes Fibroblasts Hepatocytes
# F_Cohort_1                357     829             45         205                    1424         9405          58        1879
# F_Cohort_2                358     785             24          94                     936         3666          28         492
# F_Cohort_3                296     739             13         100                     856         3375          17         265
# M_Cohort_1                209     482             46         979                     916         9560          93        2711
# M_Cohort_2                 16     620             13         104                      63         5128          80         976
# M_Cohort_3                172     322             28         100                     649         4317          17         254
# 
#             Hepatocytes_Efferocytosing HSPCs Kidney_distal_tubule Kidney_prox_tubule Lymphoid_progenitors Macrophages Mast_cells
# F_Cohort_1                       1689   267                  103                291                   85         366         50
# F_Cohort_2                       1354   210                   61                233                   79         240         56
# F_Cohort_3                        883   230                    1                  8                   45         288        322
# M_Cohort_1                       2813   231                  138                464                   77         550         98
# M_Cohort_2                       1273     0                    0                  0                    0         619        427
# M_Cohort_3                        702   226                    2                 24                   45         313        274
# 
#             Multipotent_progenitors Myeloid_progenitors Neutrophil_Progenitors Neutrophils NK_T_cells NK_T_progenitor_cells Thrombocytes
# F_Cohort_1                      60                  27                    234        1612        452                    18          158
# F_Cohort_2                      54                  27                    202        1389        215                    20           45
# F_Cohort_3                      54                  68                    233        1298        252                    16          153
# M_Cohort_1                      92                  70                    186        1229        343                    21          237
# M_Cohort_2                      40                  60                     12         200        274                    55          308
# M_Cohort_3                      40                  62                    206        1649        367                    25           56

# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(killi.atlas.sce.cl, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one sheet per subpopulation
assayNames(pb)
# [1] "B_Cell_Progenitors"         "B_cells"                    "Cholangiocytes"             "Endothelial"                "Erythrocyte_Progenitors"   
# [6] "Erythrocytes"               "Fibroblasts"                "Hepatocytes"                "Hepatocytes_Efferocytosing" "HSPCs"                     
# [11] "Kidney_distal_tubule"       "Kidney_prox_tubule"         "Lymphoid_progenitors"       "Macrophages"                "Mast_cells"                
# [16] "Multipotent_progenitors"    "Myeloid_progenitors"        "Neutrophil_Progenitors"     "Neutrophils"                "NK_T_cells"                
# [21] "NK_T_progenitor_cells"      "Thrombocytes"   

# Pseudobulk-level MDS plot
pb_mds <- pbMDS(pb)

# output global MDS
pdf(paste0(Sys.Date(),"_10x_Killi_Tissue_Atlas_Muscat_PB_MDS.pdf"))
pb_mds
dev.off()

# nb. of cells per cluster-sample
cell.per.samp.tab <- t(table(killi.atlas.sce.cl$cluster_id, killi.atlas.sce.cl$sample_id))

# cell types with at least 25 cells from every each sex/cohort sample
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab > 25, 2, sum) == 6]
# keep:
# [1] "B_cells"                    "Endothelial"                "Erythrocyte_Progenitors"    "Erythrocytes"               "Hepatocytes"               
# [6] "Hepatocytes_Efferocytosing" "Macrophages"                "Mast_cells"                 "Multipotent_progenitors"    "Myeloid_progenitors"       
# [11] "Neutrophils"                "NK_T_cells"                 "Thrombocytes"   

### extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- pb@assays@data[celltype.qc]

# some genes have "gene-" in their name: remove for processing
for (i in 1:length(counts.pb)) {
  rownames(counts.pb[[i]]) <- gsub("gene-", "",rownames(counts.pb[[i]]))
}


###############################################################################################
# 2. Use RUV to clean up batch/tissue effects on expression

# will run RUV with k = 2
#     - batch + different tissue contribution to the pseudo bulk => 2 sources
# run for the cell types with at least 25 cells from every each sex/cohort 

# Make the group object for RUV (will be identical for all of them)
differences <- makeGroups(ifelse(grepl("F_",colnames(counts.pb[[1]])),"F","M"))

# Create list object to receive clean RUV counts
ruv.cleaned.counts <- vector(mode = "list", length = length(counts.pb))
names(ruv.cleaned.counts) <- names(counts.pb)

##################################################################
######################## Try alternate pipeline: edgeR-LRT / 2023-06-08

# Create list object to receive VST normalized counts
tmm.counts <- vector(mode = "list", length = length(counts.pb))
names(tmm.counts) <- names(counts.pb)

# Create list object to receive edgeR-LRT results
edger.res.list <- vector(mode = "list", length = length(counts.pb))
names(edger.res.list) <- names(counts.pb)

# Create list object to receive edgeR-LRT results (just genes)
edger.res.list.genes <- vector(mode = "list", length = length(counts.pb))
names(edger.res.list.genes) <- names(counts.pb)


# loop over pseudobulk data
for  (i in 1:length(counts.pb)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_edgeR_LRT_Pseudobulk_",names(counts.pb)[[i]])
  
  # Do RUV, then remove genes not expressed robustly (in at least 2 samples)
  my.tmp.ruv <- (RUVs(counts.pb[[i]], rownames(counts.pb[[i]]), k=2, differences))$normalizedCounts # k = 2 RUV: batch + different tissue contribution
  my.good    <- apply(my.tmp.ruv>0,1,sum) >= 2
  
  # keep only robustly expressed genes
  ruv.cleaned.counts[[i]] <- my.tmp.ruv[my.good,]
  
  my.cols <- c(rep("deeppink",3),rep("deepskyblue",3))
  
  my.sex <- factor(ifelse(grepl("F_",colnames(counts.pb[[i]])),"F","M"), levels = c("M","F")) # so that females have positive logFC
  
  # create DGElist for edger
  edger.obj <- DGEList(counts = ruv.cleaned.counts[[i]], 
                       group = my.sex)
  
  # filter on expression
  keep      <- filterByExpr(edger.obj, min.counts = 3, min.prop = 0.5) # 50% of samples at least 1 count
  edger.obj <- edger.obj[keep,,keep.lib.sizes=FALSE]
  
  # normalize library sizes
  edger.obj <- calcNormFactors(edger.obj, method="TMM")
  design    <- model.matrix(~ my.sex)
  edger.obj <- estimateDisp(edger.obj, design)
  
  # get edger TMM normalized expression value
  tmm.counts[[i]] <- cpm(edger.obj)

  
  # MDS analysis
  mds.result <- cmdscale(1-cor(tmm.counts[[i]], method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  pdf(paste0(my.outprefix,"_TMM_MDS_plot.pdf"))
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main= paste0(names(counts.pb)[[i]]," MDS"),
       cex=3, col= my.cols, pch= 16,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  dev.off()
  
  # Perform likelihood ratio tests
  edger.fit <- glmFit(edger.obj, design)
  edger.lrt <- glmLRT(edger.fit, coef=2)
  
  # extract gene significance by edger-LRT
  res.sex <- topTags(edger.lrt, n = Inf)$table # FC in females over Males
  
  # restrict analysis to only genes, not TEs
  res.sex.genes <- res.sex[!grepl("NotFur", rownames(res.sex)), ]
  
  # store results
  edger.res.list[[i]]       <- data.frame(res.sex)
  edger.res.list.genes[[i]] <- data.frame(res.sex.genes)
  
  ### get sex dimorphic changes at FDR5
  genes.sex  <- rownames(res.sex.genes)[res.sex.genes$FDR < 0.05]
  my.num.sex <- length(genes.sex)
  
  if (my.num.sex > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste(my.outprefix,"SEX_DIM_Heatmap_FDR5_GENES.pdf", sep = "_")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb)[[i]], " sex significant (FDR<5%), ", my.num.sex, " genes")
    pheatmap::pheatmap(tmm.counts[[i]][genes.sex,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = F, scale="row",
                       main = my.heatmap.title, 
                       cellwidth = 20, 
                       border    = NA,
                       cellheight = 1 )
    dev.off()
  }
  
  # output result tables of combined analysis to text files
  my.out.ct.mat <- paste0(my.outprefix,"_TMM_log2_counts_matrix.txt")
  write.table(tmm.counts[[i]], file = my.out.ct.mat , sep = "\t" , row.names = T, quote = F)
  
  my.out.stats.sex <- paste0(my.outprefix,"_SEX_DIM_all_genes_TEs_statistics.txt")
  write.table(edger.res.list[[i]], file = my.out.stats.sex , sep = "\t" , row.names = T, quote = F)
  
  my.out.fdr5.sex <- paste0(my.outprefix,"_SEX_DIM_FDR5_genes_TEs_statistics.txt")
  write.table(edger.res.list[[i]][genes.sex,], file = my.out.fdr5.sex, sep = "\t" , row.names = T, quote = F)
  
  
  #### now do TEs
  te.sex <- rownames(res.sex)[res.sex$FDR < 0.05][grep("NotFur1-",rownames(res.sex)[res.sex$FDR < 0.05])]
  
  if (length(te.sex) > 2) {
    # heatmap drawing - only if there is at least 2 TEs
    my.heatmap.out <- paste(my.outprefix,"SEX_DIM_Heatmap_FDR5_TEs.pdf", sep = "_")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 15)
    my.heatmap.title <- paste0(names(counts.pb)[[i]], " sex significant (FDR<5%), ", length(te.sex), " TEs")
    pheatmap::pheatmap(tmm.counts[[i]][te.sex,],
                       cluster_cols = F,
                       cluster_rows = T,
                       colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
                       show_rownames = T, scale="row",
                       main = my.heatmap.title, 
                       cellwidth = 20, 
                       border = NA,
                       cellheight = 10 )
    dev.off()
  }
  
}

# save R object with all EDGER results
my.rdata.sex <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_SEX_DIM_edger_LRT_objects.RData")
save(edger.res.list, edger.res.list.genes, file = my.rdata.sex)

my.tmm.sex <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_SEX_DIM_TMM_data_objects.RData")
save(tmm.counts, file = my.tmm.sex)


get_sig_num <- function(edger.df) {
  genes.sex <- rownames(edger.df)[edger.df$FDR < 0.05]
  my.num.sex <- length(genes.sex)
  my.num.sex
}

unlist(lapply(edger.res.list, get_sig_num))
# B_cells                Endothelial    Erythrocyte_Progenitors               Erythrocytes 
# 158                         24                          8                         35 
# Hepatocytes Hepatocytes_Efferocytosing                Macrophages                 Mast_cells 
# 535                         44                        191                         31 
# Multipotent_progenitors        Myeloid_progenitors                Neutrophils                 NK_T_cells 
# 1                        122                         55                        146 
# Thrombocytes 
# 57 

unlist(lapply(edger.res.list.genes, get_sig_num))
# B_cells                Endothelial    Erythrocyte_Progenitors               Erythrocytes 
# 155                         20                          7                         29 
# Hepatocytes Hepatocytes_Efferocytosing                Macrophages                 Mast_cells 
# 519                         34                        183                         31 
# Multipotent_progenitors        Myeloid_progenitors                Neutrophils                 NK_T_cells 
# 1                        104                         54                        118 
# Thrombocytes 
# 54 

pdf(paste0(Sys.Date(),"_barplot_number_of_significant_sex_biased_genes_per_cell_Type_FDR5.pdf"), height = 5, width = 7)
par(oma=c(0.5,7,0.5,0.5))
barplot(unlist(lapply(edger.res.list.genes, get_sig_num)), las = 1, horiz = T, 
        xlab = "significant genes", col = "black", main = "Sex-biased genes (FDR < 0.05)",
        xlim = c(0,600))
box()
dev.off()

# Split by direction
get_sig_num_sex <- function(edger.df) {
  genes.sex.F <- rownames(edger.df)[bitAnd(edger.df$FDR < 0.05, edger.df$logFC >0)>0]
  genes.sex.M <- rownames(edger.df)[bitAnd(edger.df$FDR < 0.05, edger.df$logFC <0)>0]
  c(length(genes.sex.F),length(genes.sex.M))
}

my.res.df <- data.frame(matrix(0,length(edger.res.list.genes),2))
rownames(my.res.df) <- names(edger.res.list.genes)
colnames(my.res.df) <- c("F_biased","M_biased")

for (i in 1:length(edger.res.list.genes)) {
  my.res.df[i,] <- get_sig_num_sex(edger.res.list.genes[[i]])
  
}

pdf(paste0(Sys.Date(),"_barplot_number_of_significant_sex_biased_genes_per_cell_Type_FDR5_BY_SEX.pdf"), height = 5, width = 6)
par(oma=c(0.5,7,0.5,0.5))
barplot(t(as.matrix(my.res.df)), las = 1, horiz = T, beside = T, col = c("deeppink","deepskyblue"),
        xlab = "significant genes", main = "Sex-biased genes (FDR < 0.05)",
        xlim = c(0,400))
box()
dev.off()


### Try jitter plot 
# load("2023-06-08_pseudobulk_killi_cell_types_SEX_DIM_edger_LRT_objects.RData")

# code adapted from https://github.com/brunetlab/Leeman_et_al_2017/blob/master/kallisto_deseq2/Fig4A_stripplot_cell_type_colors.R
# edger.res.list

# get color palette used for UMAP
colors.cells <- read.table('2023-04-12_color_palette_annotation.txt', header = T)

# get correspondence cell type/color
ct.unique <- data.frame("label" = colors.cells$cell_type, "cols" = colors.cells$col)

# filter to EdgeR tested cell types only
ct.unique           <- ct.unique[ct.unique$label %in% names(edger.res.list),]
rownames(ct.unique) <- ct.unique$label

## Order by pvalue:
sex.results <- lapply(edger.res.list,function(x) {x[order(x$FDR),]})
n           <- sapply(sex.results, nrow)
names(n)    <- names(sex.results)

cols <- list()
xlab <- character(length = length(sex.results))
for(i in seq(along = sex.results)){
  cols[[i]] <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i <- sex.results[[i]]$FDR < 0.05
  cols[[i]][ind.sig.i] <- ct.unique[names(sex.results)[i], "cols"]
  xlab[i] <- paste(names(sex.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(sex.results)

pdf(paste0(Sys.Date(),"_stripplot_edgeR_with_cell_type_colors.pdf"), width = 6, height = 5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 13.5),
     ylim = c(-15, 20),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change (F / M)"
)
abline(h = 0)
abline(h = seq(-20, 20, by = 5)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(sex.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(sex.results[[i]])), amount = 0.2),
         y = rev(sex.results[[i]]$logFC),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]),
         cex = 0.7)
}
axis(1,
     at = 1:13,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-15, 20, 5))
box()
dev.off()
###############################################################################################



#######################
sink(file = paste(Sys.Date(),"_scRNAseq_PseudoBulk_analysis_session_Info_EDGER.txt", sep =""))
sessionInfo()
sink()

