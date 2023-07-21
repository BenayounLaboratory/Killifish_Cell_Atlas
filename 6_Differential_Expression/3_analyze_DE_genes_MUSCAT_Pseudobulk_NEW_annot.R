setwd('/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2022-02-01_Killifish_scRNAseq_ATLAS_analyses/Differential_Expression/')
options(stringsAsFactors = F)

### R version 3.6.3 (2020-02-29)

library('Seurat')         # Seurat_3.2.2 
library('muscat')         # muscat_1.0.1    
library(sctransform)      # sctransform_0.3.2 
library("singleCellTK")   # singleCellTK_1.6.0 
library('DESeq2')         # DESeq2_1.26.0  
library('RUVSeq')         # RUVSeq_1.20.0  
library(clusterProfiler)  # clusterProfiler_3.14.3
library(org.Dr.eg.db)     # org.Dr.eg.db_3.10.0
require(DOSE)             # DOSE_3.12.0     
library(ggplot2)          # ggplot2_3.3.5    
library(scales)           # scales_1.1.1 

library("ComplexHeatmap") # 
library("bitops")         # 
library(circlize)         #
library(Polychrome)

theme_set(theme_bw())   


# 2023-02-10
# Run on clean new annotation

###############################################################################################
# 0. preprocess Seurat object for use with muscat

# Load annotated cell data
load('../Seurat_merged_analysis/2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData')
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
killi.atlas.sce.cl <- prepSCE(killi.atlas.sce, 
                              cluster_id = "Annotation_v1",  # population assignments
                              group_id   = "Sex"          ,  # group IDs (ctrl/stim)
                              sample_id  = "Fish"        ,  # sample IDs (ctrl/stim.1234)
                              drop       = TRUE)             # drop all other colData columns


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

# Create list object to receive VST normalized counts
vst.counts <- vector(mode = "list", length = length(counts.pb))
names(vst.counts) <- names(counts.pb)

# Create list object to receive DESeq2 results
deseq.res.list <- vector(mode = "list", length = length(counts.pb))
names(deseq.res.list) <- names(counts.pb)

# Create list object to receive DESeq2 results (just genes)
deseq.res.list.genes <- vector(mode = "list", length = length(counts.pb))
names(deseq.res.list.genes) <- names(counts.pb)


# loop over pseudobulk data
for  (i in 1:length(counts.pb)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_DEseq2_Pseudobulk_",names(counts.pb)[[i]])
  
  # Do RUV, then remove genes not expressed robustly (in at least 2 samples)
  my.tmp.ruv <- (RUVs(counts.pb[[i]], rownames(counts.pb[[i]]), k=2, differences))$normalizedCounts # k = 2 RUV: batch + different tissue contribution
  my.good    <- apply(my.tmp.ruv>0,1,sum) >= 2
  
  # keep only robustly expressed genes
  ruv.cleaned.counts[[i]] <- my.tmp.ruv[my.good,]
  
  my.cols <- c(rep("deeppink",3),rep("deepskyblue",3))
  
  # build design matrix
  dataDesign = data.frame(row.names = colnames( ruv.cleaned.counts[[i]] ), 
                          sex = ifelse(grepl("F_",colnames(counts.pb[[i]])),"F","M"))
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = ruv.cleaned.counts[[i]],
                                colData = dataDesign,
                                design = ~ sex)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  # get DESeq2 normalized expression value
  vst.counts[[i]] <- getVarianceStabilizedData(dds.deseq)
  
  # MDS analysis
  mds.result <- cmdscale(1-cor(vst.counts[[i]],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  pdf(paste0(my.outprefix,"_MDS_plot.pdf"))
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main= paste0(names(counts.pb)[[i]]," MDS"),
       cex=3, col= my.cols, pch= 16,
       cex.lab = 1.25,
       cex.axis = 1.25, las = 1)
  dev.off()
  
  # extract gene significance by DEseq2
  res.sex <- results(dds.deseq, contrast = c("sex","F","M")) # FC in females over Males
  
  # exclude genes with NA FDR value
  res.sex <- res.sex[!is.na(res.sex$padj),]
  
  # restrict analysis to only genes, not TEs
  res.sex.genes <- res.sex[!grepl("NotFur", rownames(res.sex)), ]
  
  # store results
  deseq.res.list[[i]] <- data.frame(res.sex)
  deseq.res.list.genes[[i]] <- data.frame(res.sex.genes)
  
  ### get sex dimorphic changes at FDR5
  genes.sex <- rownames(res.sex.genes)[res.sex.genes$padj < 0.05]
  my.num.sex <- length(genes.sex)
  
  if (my.num.sex > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste(my.outprefix,"SEX_DIM_Heatmap_FDR5_GENES.pdf", sep = "_")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 10)
    my.heatmap.title <- paste0(names(counts.pb)[[i]], " sex significant (FDR<5%), ", my.num.sex, " genes")
    pheatmap::pheatmap(vst.counts[[i]][genes.sex,],
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
  my.out.ct.mat <- paste0(my.outprefix,"_VST_log2_counts_matrix.txt")
  write.table(vst.counts[[i]], file = my.out.ct.mat , sep = "\t" , row.names = T, quote = F)
  
  my.out.stats.sex <- paste0(my.outprefix,"_SEX_DIM_all_genes_TEs_statistics.txt")
  write.table(deseq.res.list[[i]], file = my.out.stats.sex , sep = "\t" , row.names = T, quote = F)
  
  my.out.fdr5.sex <- paste0(my.outprefix,"_SEX_DIM_FDR5_genes_TEs_statistics.txt")
  write.table(deseq.res.list[[i]][genes.sex,], file = my.out.fdr5.sex, sep = "\t" , row.names = T, quote = F)
  
  
  #### now do TEs
  te.sex <- rownames(res.sex)[res.sex$padj < 0.05][grep("NotFur1-",rownames(res.sex)[res.sex$padj < 0.05])]
  
  if (length(te.sex) > 2) {
    # heatmap drawing - only if there is at least 2 gene
    my.heatmap.out <- paste(my.outprefix,"SEX_DIM_Heatmap_FDR5_TEs.pdf", sep = "_")
    
    pdf(my.heatmap.out, onefile = F, height = 10, width = 15)
    my.heatmap.title <- paste0(names(counts.pb)[[i]], " sex significant (FDR<5%), ", length(te.sex), " TEs")
    pheatmap::pheatmap(vst.counts[[i]][te.sex,],
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

# save R object with all DEseq2 results
my.rdata.sex <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_SEX_DIM_DEseq2_objects.RData")
save(deseq.res.list, deseq.res.list.genes, file = my.rdata.sex)

my.vst.sex <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_SEX_DIM_VST_data_objects.RData")
save(vst.counts, file = my.vst.sex)


get_sig_num <- function(deseq2.df) {
  genes.sex <- rownames(deseq2.df)[deseq2.df$padj < 0.05]
  my.num.sex <- length(genes.sex)
  my.num.sex
}

unlist(lapply(deseq.res.list, get_sig_num))
# B_cells                Endothelial    Erythrocyte_Progenitors               Erythrocytes                Hepatocytes 
# 81                         11                          1                         15                        311 
# Hepatocytes_Efferocytosing                Macrophages                 Mast_cells    Multipotent_progenitors        Myeloid_progenitors 
# 17                        166                         26                          0                         20 
# Neutrophils                 NK_T_cells               Thrombocytes 
# 34                         30                         31 


unlist(lapply(deseq.res.list.genes, get_sig_num))
# B_cells                Endothelial    Erythrocyte_Progenitors               Erythrocytes                Hepatocytes 
# 81                         11                          1                         13                        306 
# Hepatocytes_Efferocytosing                Macrophages                 Mast_cells    Multipotent_progenitors        Myeloid_progenitors 
# 17                        162                         26                          0                         19 
# Neutrophils                 NK_T_cells               Thrombocytes 
# 34                         24                         31 

pdf(paste0(Sys.Date(),"_barplot_number_of_significant_sex_biased_genes_per_cell_Type_FDR5.pdf"), height = 5, width = 7)
par(oma=c(0.5,7,0.5,0.5))
barplot(unlist(lapply(deseq.res.list.genes, get_sig_num)), las = 1, horiz = T, 
        xlab = "significant genes", col = "black", main = "Sex-biased genes (FDR < 0.05)",
        xlim = c(0,400))
box()
dev.off()

# Split by direction
get_sig_num_sex <- function(deseq2.df) {
  genes.sex.F <- rownames(deseq2.df)[bitAnd(deseq2.df$padj < 0.05, deseq2.df$log2FoldChange >0)>0]
  genes.sex.M <- rownames(deseq2.df)[bitAnd(deseq2.df$padj < 0.05, deseq2.df$log2FoldChange <0)>0]
  c(length(genes.sex.F),length(genes.sex.M))
}

my.res.df <- data.frame(matrix(0,length(deseq.res.list.genes),2))
rownames(my.res.df) <- names(deseq.res.list.genes)
colnames(my.res.df) <- c("F_biased","M_biased")

for (i in 1:length(deseq.res.list.genes)) {
  my.res.df[i,] <- get_sig_num_sex(deseq.res.list.genes[[i]])
  
}

pdf(paste0(Sys.Date(),"_barplot_number_of_significant_sex_biased_genes_per_cell_Type_FDR5_BY_SEX.pdf"), height = 5, width = 6)
par(oma=c(0.5,7,0.5,0.5))
barplot(t(as.matrix(my.res.df)), las = 1, horiz = T, beside = T, col = c("deeppink","deepskyblue"),
        xlab = "significant genes", main = "Sex-biased genes (FDR < 0.05)",
        xlim = c(0,200))
box()
dev.off()


##########################################################################################
### Try jitter plot (2023-04-10)
load("2023-02-10_pseudobulk_killi_cell_types_SEX_DIM_DEseq2_objects.RData")

# code adapted from https://github.com/brunetlab/Leeman_et_al_2017/blob/master/kallisto_deseq2/Fig4A_stripplot_cell_type_colors.R
# deseq.res.list

# get color palette used for UMAP
colors.cells <- read.table('../Seurat_merged_analysis/2023-04-12_color_palette_annotation.txt', header = T)

# get correspondence cell type/color
ct.unique <- data.frame("label" = colors.cells$cell_type, "cols" = colors.cells$col)

# filter to DEseq2 tested cell types only
ct.unique           <- ct.unique[ct.unique$label %in% names(deseq.res.list),]
rownames(ct.unique) <- ct.unique$label

## Order by pvalue:
sex.results <- lapply(deseq.res.list,function(x) {x[order(x$padj),]})
n        <- sapply(sex.results, nrow)
names(n) <- names(sex.results)

cols <- list()
xlab <- character(length = length(sex.results))
for(i in seq(along = sex.results)){
  cols[[i]] <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i <- sex.results[[i]]$padj < 0.05
  cols[[i]][ind.sig.i] <- ct.unique[names(sex.results)[i], "cols"]
  xlab[i] <- paste(names(sex.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(sex.results)

pdf(paste0(Sys.Date(),"_stripplot_DESeq2_with_cell_type_colors.pdf"), width = 6, height = 5)
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
         y = rev(sex.results[[i]]$log2FoldChange),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]))
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


###############################################################################################
# 3. clusterprofiler GO BP analysis

# Create list object to receive clean RUV counts
go.results <- vector(mode = "list", length = length(deseq.res.list.genes))
names(go.results) <- names(deseq.res.list.genes)

# Read in BLAST homology file for killifish/zebrafish (best zebrafish hit to killifish to get conversion)
zebra.homol <- read.csv("../Zebrafish_alignment/2022-03-15_Zebrafish_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-5_Minimal_HOMOLOGY_TABLE_REV.txt", sep = "\t", header = T)

# Loop over DEseq2 results
for (i in 1:length(deseq.res.list.genes)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_ClusterProfiler_GSEA_GOBP_",names(deseq.res.list.genes)[[i]])
  
  # Nfur Data
  sc.res             <- deseq.res.list.genes[[i]]
  sc.res$Nfur_Symbol <- rownames(sc.res)
  
  # Get Zebrafish homolog information symbols based on the BLAST results file using org.Dr.eg.db package
  # Some Ids will fail to map and will be ignored
  dataZ     = merge(unique(zebra.homol[,c("Nfur_Symbol","DanRer_Symbol")]), sc.res, by = "Nfur_Symbol") 
  head(dataZ)
  
  # There can be duplicate values because of paralogs, I take average of those for log2FoldChange
  data.unique = aggregate(dataZ[,"log2FoldChange"], list(dataZ$DanRer_Symbol), mean)
  colnames(data.unique) = c("zebrafish", "log2FoldChange")
  data.unique <- data.unique[data.unique != "",] # discard if no zebrafish homolog
  head(data.unique)
  
  
  # generate and sort the gene list based on log2FoldChange in decreasing order. 
  geneList = data.unique[,2]  # gene list for GO 
  names(geneList) = as.character(data.unique[,1]) # with zebrafish gene symbols as names
  geneList = sort(geneList, decreasing = TRUE)
  
  
  # Gene Ontology ------------------------------------------------------------------------------------------------------------------------------------
  go.bp.gsea <- gseGO(geneList     = geneList,
                      OrgDb        = org.Dr.eg.db,
                      keyType      = 'SYMBOL',
                      ont          = "BP",
                      nPerm        = 1000,
                      minGSSize    = 25,
                      maxGSSize    = 5000)
  View(go.bp.gsea@result)
  
  # write results to file
  write.table(go.bp.gsea@result, file = paste0(my.outprefix,"_GSEA_Analysis_FDR5.txt"), quote = F, sep = "\t")
  
  # save results into object
  go.results[[i]] <- go.bp.gsea@result
  
}

# save R object with all DEseq2 results
my.rdata.GO.sex <- paste0(Sys.Date(),"_pseudobulk_killi_cell_types_SEX_DIM_GOBP_GSEA.RData")
save(go.results, file = my.rdata.GO.sex)


### Make bubble plots
source('2022-08-19_Plot_bubble_chart_function.R')

for(i in 1:length(go.results)) {
  if(nrow(go.results[[i]])>0) { # only if something passed significance
    get_sex_bubble_plot(names(go.results)[i], go.results[[i]])
  }
}


###################### Summary plot across cell types ######################  

# get data from significant FDR 0.05
my.tissues.go <- vector(length=length(go.results), mode="list")
names(my.tissues.go) <- names(go.results)
my.pathways <- c()

for ( i in 1:length(go.results)) {
  my.GO.terms <- paste(go.results[[i]]$ID, go.results[[i]]$Description)
  go.results[[i]]$GO_Term <- my.GO.terms
  my.pathways <- unique(c(my.pathways,my.GO.terms))
}


############
# prepapre output data
# p-val matrix
my.matrix <- matrix(0,length(my.pathways),length(go.results)) # default: -log10(1) pval == 0 no enrichment

# Enrichment matrix
my.matrix2 <- matrix(0,length(my.pathways),length(go.results)) # initialize with Enrichment = 0 if no enrich

# matrix with record of significance
my.matrix3 <- matrix(0,length(my.pathways),length(go.results)) # to get sigificant pathways

colnames(my.matrix)  <- names(go.results)
colnames(my.matrix2) <- names(go.results)
colnames(my.matrix3) <- names(go.results)
rownames(my.matrix)  <- my.pathways
rownames(my.matrix2) <- my.pathways
rownames(my.matrix3) <- my.pathways

# collect data from clusterprofiler run
for (i in 1:length(my.pathways)) {
  for (j in 1:length(go.results)) { # tissues 
    # determine position
    my.id <- which(go.results[[j]]$GO_Term %in% my.pathways[i])
    
    if(length(my.id) == 1) { # if was present in this tissue
      # extract stats: FDR, NES and presence
      my.matrix[i,j] <- -log10(go.results[[j]]$p.adjust[my.id]) # log(0) is undefined
      my.matrix2[i,j] <- go.results[[j]]$NES[my.id]
      my.matrix3[i,j] <- 1
    }
  }
}


# find pathways significant in 8 or more cell types [13 passed threshold]
my.sigs         <- apply(my.matrix3,1,sum) >= 6    
sum(my.sigs) # 15
my.res.enrich   <- data.frame(my.matrix2[my.sigs,])
my.pval.enrich  <- data.frame(my.matrix[my.sigs,])

# sort by average change
my.average     <- apply(my.res.enrich,1,mean)
my.sorted      <- sort(my.average,index.return=T,decreasing=T)
my.res.enrich2 <- my.res.enrich[my.sorted$ix,]

my.pval.enrich2 <- data.frame(my.pval.enrich[my.sorted$ix,])

my.txtname <- paste(Sys.Date(),"Enrichment_table_scRNAseq_killi_SexDim_GOBP_pathways_significant_in_8_or_more.txt", sep="_")
write.table(my.res.enrich2,file=my.txtname,sep="\t",quote=F)

my.res.enrich2$Pathnames <- rownames(my.res.enrich2)

# format for ggplot
my.res.enrich3 <- cbind(my.res.enrich2[,c('Pathnames',names(go.results)[1])],rep(names(go.results)[1],dim(my.res.enrich2)[1]),my.pval.enrich2[,names(go.results)[1]])
colnames(my.res.enrich3) <- c('Pathnames','NES','condition','minusLog10Pval')
for ( h in 2:length(names(go.results))) {
  my.new <- cbind(my.res.enrich2[,c('Pathnames',names(go.results)[h])],rep(names(go.results)[h],dim(my.res.enrich2)[1]),my.pval.enrich2[,names(go.results)[h]])
  colnames(my.new) <- colnames(my.res.enrich3)
  my.res.enrich3 <- rbind(my.res.enrich3, 
                          my.new)
  
}

# Female/Male color scale
my.max <-  3 # add a value in case there are not terms biased in one sex
my.min <- -3 # add a value in case there are not terms biased in one sex
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector.sex <- c("deepskyblue","lightskyblue","lightskyblue1","lightcyan","white","lavenderblush","plum1","orchid1","deeppink")

# to preserve the wanted order
my.res.enrich3$condition <- factor(my.res.enrich3$condition, levels = unique(my.res.enrich3$condition))
my.res.enrich3$Pathnames <- factor(my.res.enrich3$Pathnames, levels = rev(unique(my.res.enrich3$Pathnames)))

my.pdfname <- paste(Sys.Date(),"Enrichment_BALLOON_plot_scRNAseq_killi_SexDim_GOBP_pathways_significant_in_6_or_more.pdf", sep="_")

pdf(my.pdfname, onefile=F, height = 5, width=8)
my.plot <- ggplot(my.res.enrich3,aes(x=condition,y=Pathnames,colour=NES,size=minusLog10Pval))+ theme_bw()+ geom_point(shape = 16) 
my.plot <- my.plot + ggtitle("Sex biased pathways") + labs(x = "Tissue/condition", y = "Gene Set")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.sex,na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(-3, 3) )
my.plot <- my.plot + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(my.plot)
dev.off()  
#################################################################################################################################################################



###############################################################################################
# 4. cross tissue comparison of genes with sex-dimorphic gene expression

# get the common expressed background for comparison
my.exp.genes <- rownames(killi.atlas.sce.cl)

### 
my.logFC.data <- data.frame(matrix(NA,length(my.exp.genes),length(deseq.res.list.genes)))
colnames(my.logFC.data) <- names(deseq.res.list.genes)
rownames(my.logFC.data) <- my.exp.genes

my.FDR.data <- data.frame(matrix(NA,length(my.exp.genes),length(deseq.res.list.genes)))
colnames(my.FDR.data) <- names(deseq.res.list.genes)
rownames(my.FDR.data) <- my.exp.genes

for (i in 1:length(deseq.res.list.genes)) {
  for (j in 1:length(my.exp.genes)) {
    my.idx <- rownames(deseq.res.list.genes[[i]]) %in% my.exp.genes[j]
    
    if (sum(my.idx) > 0) {
      
      my.logFC.data[j,i] <- deseq.res.list.genes[[i]]$log2FoldChange[my.idx]
      my.FDR.data[j,i]   <- deseq.res.list.genes[[i]]$padj[my.idx]
      
    }
  }
}

# get spearman rank correlation
my.cors <-cor(my.logFC.data, method = 'spearman', use = "complete.obs")

pdf(paste0(Sys.Date(),"_correlation_of_FC_heatmap.pdf"))
Heatmap(my.cors, border = T, rect_gp = gpar(col = "grey", lwd = 0.5),
        column_title = "log2FC correlation across cell types in Female vs. Males")
dev.off()


#########
source('2022-08-19_Plot_Jaccard_function.R')
##
jacc.5 <- calc_jaccard(my.logFC.data, my.FDR.data, 0.05)

pdf(paste0(Sys.Date(),"_Jaccard_Index_FDR5.pdf"))
Heatmap(jacc.5[[1]], 
        col = colorRamp2(c(0,1), c("white","deeppink"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5)  ,
        column_title = "Jaccard Female-biased")
Heatmap(jacc.5[[2]], 
        col = colorRamp2(c(0,1), c("white","deepskyblue"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        column_title = "Jaccard Male-biased")
dev.off()

##
jacc.10 <- calc_jaccard(my.logFC.data, my.FDR.data, 0.1)

pdf(paste0(Sys.Date(),"_Jaccard_Index_FDR10.pdf"))
Heatmap(jacc.10[[1]], 
        col = colorRamp2(c(0,1), c("white","deeppink"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5)  ,
        column_title = "Jaccard Female-biased")
Heatmap(jacc.10[[2]], 
        col = colorRamp2(c(0,1), c("white","deepskyblue"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        column_title = "Jaccard Male-biased")
dev.off()



# genes reccurently female or male biased
my.up  <- my.logFC.data > 0
my.sig <- my.FDR.data   < 0.05

up.genes   <- data.frame(matrix(0,length(my.exp.genes),length(deseq.res.list.genes)))
colnames(up.genes) <- names(deseq.res.list.genes)
rownames(up.genes) <- my.exp.genes

down.genes <- data.frame(matrix(0,length(my.exp.genes),length(deseq.res.list.genes)))
colnames(down.genes) <- names(deseq.res.list.genes)
rownames(down.genes) <- my.exp.genes

for (i in 1:length(deseq.res.list.genes)) {
  for (j in 1:length(my.exp.genes)) {
    
    up.genes[j,i]   <- bitAnd(my.up[j,i], my.sig[j,i])
    down.genes[j,i] <- bitAnd(!my.up[j,i], my.sig[j,i])
    
  }
}
up.genes[is.na(up.genes)]     <- 0
down.genes[is.na(down.genes)] <- 0

# get top 3 most recurrently sex biased

# rownames(up.genes)[apply(up.genes,1,sum)>8]
# # "LOC107388898" zona pellucida sperm-binding protein 3-like
# # "LOC107388899" zona pellucida sperm-binding protein 4-like
# # "LOC107392619" uncharacterized LOC107392619
# # "LOC107392618" vitellogenin-1-like
# # "LOC107372943" zona pellucida sperm-binding protein 4-like
# # "LOC107373895" LOC107373895 heat shock 70 kDa protein 1 
# # "LOC107373896" uncharacterized LOC107373896 
# 
# rownames(down.genes)[apply(down.genes,1,sum)>4]
# # "hpx"

sort(apply(up.genes,1,sum), decreasing = T)[1:3]
# LOC107388898 LOC107392619 LOC107373896 
# 10            9            9 
# # "LOC107388898" zona pellucida sperm-binding protein 3-like
# # "LOC107392619" uncharacterized LOC107392619
# # "LOC107373896" uncharacterized LOC107373896 

sort(apply(down.genes,1,sum), decreasing = T)[1:3]
# hpx LOC107395221 LOC107375611 
# 5            4            3 
# LOC107395221 apoa1a apolipoprotein A-Ia [ Nothobranchius furzeri (turquoise killifish) ]
# LOC107375611 aldob aldolase b, fructose-bisphosphate [ Nothobranchius furzeri (turquoise killifish) ]

################################################################
## load vst counts for plotting
load('2023-02-10_pseudobulk_killi_cell_types_SEX_DIM_VST_data_objects.RData')

################################################################
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

################
hpx.data <- matrix(NA,length(vst.counts),6)
rownames(hpx.data) <- names(vst.counts)
colnames(hpx.data) <- colnames(vst.counts[[1]])

for (i in 1:length(vst.counts)) {
  if( sum(rownames(vst.counts[[i]]) %in% "hpx") == 1 ) {
    hpx.data[i,] <- vst.counts[[i]]["hpx",]/colSums(vst.counts[[i]])
  }
}


pdf(paste0(Sys.Date(),"_hpx_expression_VST.pdf"))
Heatmap(scale_rows(hpx.data), 
        col = colorRamp2(c(-2,0,2), c("darkblue", "white", "red"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        cluster_rows = F,
        cluster_columns = F,
        column_title = "VST hpx expression")
dev.off()


################
LOC107395221.data <- matrix(NA,length(vst.counts),6)
rownames(LOC107395221.data) <- names(vst.counts)
colnames(LOC107395221.data) <- colnames(vst.counts[[1]])

for (i in 1:length(vst.counts)) {
  if( sum(rownames(vst.counts[[i]]) %in% "LOC107395221") == 1 ) {
    LOC107395221.data[i,] <- vst.counts[[i]]["LOC107395221",]/colSums(vst.counts[[i]])
  }
}


pdf(paste0(Sys.Date(),"_LOC107395221_apoa1a_expression_VST.pdf"))
Heatmap(scale_rows(LOC107395221.data), 
        col = colorRamp2(c(-2,0,2), c("darkblue", "white", "red"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        cluster_rows = F,
        cluster_columns = F,
        column_title = "VST LOC107395221/apoa1a expression")
dev.off()

################
LOC107375611.data <- matrix(NA,length(vst.counts),6)
rownames(LOC107375611.data) <- names(vst.counts)
colnames(LOC107375611.data) <- colnames(vst.counts[[1]])

for (i in 1:length(vst.counts)) {
  if( sum(rownames(vst.counts[[i]]) %in% "LOC107375611") == 1 ) {
    LOC107375611.data[i,] <- vst.counts[[i]]["LOC107375611",]/colSums(vst.counts[[i]])
  }
}


pdf(paste0(Sys.Date(),"_LOC107375611_aldob_expression_VST.pdf"))
Heatmap(scale_rows(LOC107375611.data), 
        col = colorRamp2(c(-2,0,2), c("darkblue", "white", "red"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        cluster_rows = F,
        cluster_columns = F,
        column_title = "VST LOC107375611/aldob expression")
dev.off()





################
LOC107373896.data <- matrix(0,length(vst.counts),6)
rownames(LOC107373896.data) <- names(vst.counts)
colnames(LOC107373896.data) <- colnames(vst.counts[[1]])

for (i in 1:length(vst.counts)) {
  if( sum(rownames(vst.counts[[i]]) %in% "LOC107373896") == 1 ) {
    LOC107373896.data[i,] <- vst.counts[[i]]["LOC107373896",]/colSums(vst.counts[[i]])
  }
}


pdf(paste0(Sys.Date(),"_LOC107373896_ncfem1_expression_VST.pdf"))
Heatmap(scale_rows(LOC107373896.data), 
        col = colorRamp2(c(-2,0,2), c("darkblue", "white", "red"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        cluster_rows = F,
        cluster_columns = F,
        column_title = "VST LOC107373896/ncfem1 expression")
dev.off()


################
LOC107392619.data <- matrix(NA,length(vst.counts),6)
rownames(LOC107392619.data) <- names(vst.counts)
colnames(LOC107392619.data) <- colnames(vst.counts[[1]])

for (i in 1:length(vst.counts)) {
  if( sum(rownames(vst.counts[[i]]) %in% "LOC107392619") == 1 ) {
    LOC107392619.data[i,] <- vst.counts[[i]]["LOC107392619",]/colSums(vst.counts[[i]])
  }
}


pdf(paste0(Sys.Date(),"_LOC107392619_expression_VST.pdf"))
Heatmap(scale_rows(LOC107392619.data), 
        col = colorRamp2(c(-2,0,2), c("darkblue", "white", "red"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        cluster_rows = F,
        cluster_columns = F,
        column_title = "VST LOC107392619 expression")
dev.off()


################
LOC107388898.data <- matrix(NA,length(vst.counts),6)
rownames(LOC107388898.data) <- names(vst.counts)
colnames(LOC107388898.data) <- colnames(vst.counts[[1]])

for (i in 1:length(vst.counts)) {
  if( sum(rownames(vst.counts[[i]]) %in% "LOC107388898") == 1 ) {
    LOC107388898.data[i,] <- vst.counts[[i]]["LOC107388898",]/colSums(vst.counts[[i]])
  }
}


pdf(paste0(Sys.Date(),"_LOC107388898_expression_VST.pdf"))
Heatmap(scale_rows(LOC107388898.data), 
        col = colorRamp2(c(-2,0,2), c("darkblue", "white", "red"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        cluster_rows = F,
        cluster_columns = F,
        column_title = "VST LOC107388898/zp3 expression")
dev.off()



#################################################################################################################

#######################
sink(file = paste(Sys.Date(),"_scRNAseq_PseudoBulk_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

