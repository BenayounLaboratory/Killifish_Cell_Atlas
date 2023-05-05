setwd('/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2022-01-20_Killifish_tissue_scRNAseq_cohort3/Seurat_Doublet_Cleaning_SoupX/')
options(stringsAsFactors = F)

library('Seurat')
library(bitops)
library(sctransform)
library(SingleR)
library(clustree)
library(scales)
library(dplyr)
library(DoubletFinder)

# Add cxds_bcds_hybrid related packages
library(scds)
library(scater)
library(bitops)

# removal of ambient RNA
library('SoupX')


#####################################################################################################################
#### 1. Read and process data from cellranger for SoupX

# Calculate and clean the contribution of ambient RNA with SoupX
# https://github.com/constantAmateur/SoupX

# read 10x libraries cell ranger gene barcode matrices for soupx
Female_Kidney_c3   <- load10X('../CellRanger/Female_Kidney_3_FishTEDB_NR/outs/')
Female_Liver_c3    <- load10X('../CellRanger/Female_Liver_3_FishTEDB_NR/outs/')
Female_Spleen_c3   <- load10X('../CellRanger/Female_Spleen_3_FishTEDB_NR/outs/')
Male_Kidney_c3     <- load10X('../CellRanger/Male_Kidney_3_FishTEDB_NR/outs/')
Male_Liver_c3      <- load10X('../CellRanger/Male_Liver_3_FishTEDB_NR/outs/')
Male_Spleen_c3     <- load10X('../CellRanger/Male_Spleen_3_FishTEDB_NR/outs/')

# estimate RNA soup from bottom 10% droplets
Female_Kidney_c3     <- autoEstCont(Female_Kidney_c3    ) # Estimated global rho of 0.01
Female_Liver_c3      <- autoEstCont(Female_Liver_c3     ) # Estimated global rho of 0.01
Female_Spleen_c3     <- autoEstCont(Female_Spleen_c3    ) # Estimated global rho of 0.01
Male_Kidney_c3       <- autoEstCont(Male_Kidney_c3      ) # Estimated global rho of 0.01
Male_Liver_c3        <- autoEstCont(Male_Liver_c3       ) # Estimated global rho of 0.01
Male_Spleen_c3       <- autoEstCont(Male_Spleen_c3      ) # Estimated global rho of 0.01

# adjust counts based on RNA soup (and round)
out.Female_Kidney_c3  <- adjustCounts(Female_Kidney_c3  , roundToInt=TRUE)
out.Female_Liver_c3   <- adjustCounts(Female_Liver_c3   , roundToInt=TRUE)
out.Female_Spleen_c3  <- adjustCounts(Female_Spleen_c3  , roundToInt=TRUE)
out.Male_Kidney_c3    <- adjustCounts(Male_Kidney_c3    , roundToInt=TRUE)
out.Male_Liver_c3     <- adjustCounts(Male_Liver_c3     , roundToInt=TRUE)
out.Male_Spleen_c3    <- adjustCounts(Male_Spleen_c3    , roundToInt=TRUE)

# get seurat objects
seurat.Female_Kidney_c3   <- CreateSeuratObject( out.Female_Kidney_c3  )
seurat.Female_Liver_c3    <- CreateSeuratObject( out.Female_Liver_c3   )
seurat.Female_Spleen_c3   <- CreateSeuratObject( out.Female_Spleen_c3  )
seurat.Male_Kidney_c3     <- CreateSeuratObject( out.Male_Kidney_c3    )
seurat.Male_Liver_c3      <- CreateSeuratObject( out.Male_Liver_c3     )
seurat.Male_Spleen_c3     <- CreateSeuratObject( out.Male_Spleen_c3    )

# Merge objects for the cohort
# Merge Seurat objects
killi.combined.c3 <- merge(seurat.Female_Kidney_c3,
                       y =  c(seurat.Female_Liver_c3  ,
                              seurat.Female_Spleen_c3 ,
                              seurat.Male_Kidney_c3   ,
                              seurat.Male_Liver_c3    ,
                              seurat.Male_Spleen_c3   ),
                       add.cell.ids = c("Female_Kidney_Cohort3" ,
                                        "Female_Liver_Cohort3"  ,
                                        "Female_Spleen_Cohort3" ,
                                        "Male_Kidney_Cohort3"   ,
                                        "Male_Liver_Cohort3"    ,
                                        "Male_Spleen_Cohort3"   ),
                       project = "10x_killi_cohort3")

killi.combined.c3
# An object of class Seurat 
# 27834 features across 23798 samples within 1 assay 
# Active assay: RNA (27834 features, 0 variable features)

############################## metadata ##############################
# create condition label
my.Female_Kidney_Cohort3  <- grep("Female_Kidney_Cohort3" , colnames(killi.combined.c3@assays$RNA))
my.Female_Liver_Cohort3   <- grep("Female_Liver_Cohort3"  , colnames(killi.combined.c3@assays$RNA))
my.Female_Spleen_Cohort3  <- grep("Female_Spleen_Cohort3" , colnames(killi.combined.c3@assays$RNA))
my.Male_Kidney_Cohort3    <- grep("Male_Kidney_Cohort3"   , colnames(killi.combined.c3@assays$RNA))
my.Male_Liver_Cohort3     <- grep("Male_Liver_Cohort3"    , colnames(killi.combined.c3@assays$RNA))
my.Male_Spleen_Cohort3    <- grep("Male_Spleen_Cohort3"   , colnames(killi.combined.c3@assays$RNA))

#####
Condition <- rep("NA", length(colnames(killi.combined.c3@assays$RNA)))
Condition[ my.Female_Kidney_Cohort3  ]   <- "Female_Kidney_Cohort3" 
Condition[ my.Female_Liver_Cohort3   ]   <- "Female_Liver_Cohort3"  
Condition[ my.Female_Spleen_Cohort3  ]   <- "Female_Spleen_Cohort3" 
Condition[ my.Male_Kidney_Cohort3    ]   <- "Male_Kidney_Cohort3"   
Condition[ my.Male_Liver_Cohort3     ]   <- "Male_Liver_Cohort3"    
Condition[ my.Male_Spleen_Cohort3    ]   <- "Male_Spleen_Cohort3"   
Condition <- data.frame(Condition)
rownames(Condition) <- colnames(killi.combined.c3@assays$RNA)

#####
Sex <- rep("NA", length(colnames(killi.combined.c3@assays$RNA)))
Sex[ my.Female_Kidney_Cohort3  ]   <- "F" 
Sex[ my.Female_Liver_Cohort3   ]   <- "F"
Sex[ my.Female_Spleen_Cohort3  ]   <- "F" 
Sex[ my.Male_Kidney_Cohort3    ]   <- "M" 
Sex[ my.Male_Liver_Cohort3     ]   <- "M" 
Sex[ my.Male_Spleen_Cohort3    ]   <- "M"
Sex <- data.frame(Sex)
rownames(Sex) <- colnames(killi.combined.c3@assays$RNA)

#####
Tissue <- rep("NA", length(colnames(killi.combined.c3@assays$RNA)))
Tissue[ my.Female_Kidney_Cohort3  ]   <- "Kidney" 
Tissue[ my.Female_Liver_Cohort3   ]   <- "Liver"
Tissue[ my.Female_Spleen_Cohort3  ]   <- "Spleen" 
Tissue[ my.Male_Kidney_Cohort3    ]   <- "Kidney" 
Tissue[ my.Male_Liver_Cohort3     ]   <- "Liver"
Tissue[ my.Male_Spleen_Cohort3    ]   <- "Spleen" 
Tissue <- data.frame(Tissue)
rownames(Tissue) <- colnames(killi.combined.c3@assays$RNA)


#####
Batch <- rep("Cohort_3", length(colnames(killi.combined.c3@assays$RNA)))
Batch <- data.frame(Batch)
rownames(Batch) <- colnames(killi.combined.c3@assays$RNA)


# update Seurat
killi.combined.c3 <- AddMetaData(object = killi.combined.c3, metadata = as.vector(Condition), col.name = "Condition")
killi.combined.c3 <- AddMetaData(object = killi.combined.c3, metadata = as.vector(Sex)      , col.name = "Sex"      )
killi.combined.c3 <- AddMetaData(object = killi.combined.c3, metadata = as.vector(Tissue)   , col.name = "Tissue"   )
killi.combined.c3 <- AddMetaData(object = killi.combined.c3, metadata = as.vector(Batch)    , col.name = "Batch"   )

# https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
min.value = 0
min.cells = 50
genes.use <- rownames(killi.combined.c3@assays$RNA)
num.cells <- Matrix::rowSums(killi.combined.c3@assays$RNA@counts > min.value)
genes.use <- names(num.cells[which(num.cells >= min.cells)])

# remove low/null genes
my.killi.tissues <- subset(killi.combined.c3, features = genes.use)
my.killi.tissues
# An object of class Seurat 
# 15780 features across 23798 samples within 1 assay 
# Active assay: RNA (15780 features, 0 variable features)

################################################################################################################################################################
#### 2. QC on mitochondrial reads
# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
my.killi.tissues[["percent.mito"]] <- PercentageFeatureSet(my.killi.tissues, pattern = "^MT-")
head(my.killi.tissues@meta.data)

pdf(paste(Sys.Date(),"Killifish_Tissue_Cohort3violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = my.killi.tissues, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(my.killi.tissues, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(my.killi.tissues, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"Killifish_Tissue_Cohort3QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2))
dev.off()

# filter dead/low Q cells
my.killi.tissues <- subset(my.killi.tissues, subset = nFeature_RNA > 250 & percent.mito < 30 & nFeature_RNA < 7000)
my.killi.tissues
# An object of class Seurat 
# 15780 features across 21747 samples within 1 assay 
# Active assay: RNA (15780 features, 0 variable features)

################################################################################################################################################################
#### 3. Store info on biological origin of sample

head(my.killi.tissues@meta.data)
# orig.ident nCount_RNA nFeature_RNA             Condition Sex Tissue    Batch percent.mito
# Female_Kidney_Cohort3_AAACCCAAGCTGTCCG-1 SeuratProject       6137         1333 Female_Kidney_Cohort3   F Kidney Cohort_3    1.6946391
# Female_Kidney_Cohort3_AAACCCACAAAGGCTG-1 SeuratProject       6197          770 Female_Kidney_Cohort3   F Kidney Cohort_3    2.9207681
# Female_Kidney_Cohort3_AAACCCACATGAAGCG-1 SeuratProject       4356          812 Female_Kidney_Cohort3   F Kidney Cohort_3    2.7318641
# Female_Kidney_Cohort3_AAACCCAGTCAAGGCA-1 SeuratProject      23955         3167 Female_Kidney_Cohort3   F Kidney Cohort_3    1.4068044
# Female_Kidney_Cohort3_AAACCCAGTTGTCTAG-1 SeuratProject      10523          642 Female_Kidney_Cohort3   F Kidney Cohort_3    0.6462036
# Female_Kidney_Cohort3_AAACCCATCATGCCGG-1 SeuratProject      24954         1557 Female_Kidney_Cohort3   F Kidney Cohort_3    1.0138655

table(my.killi.tissues@meta.data$Condition)
# Female_Kidney_Cohort3  Female_Liver_Cohort3 Female_Spleen_Cohort3   Male_Kidney_Cohort3    Male_Liver_Cohort3   Male_Spleen_Cohort3 
#                   4650                  2786                  3215                  4405                  2773                  3918 

################################################################################################################################################################
#### 4. Normalizing the data
# global-scaling normalization method ???LogNormalize??? that normalizes the gene expression measurements for each cell by the total expression,
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
my.killi.tissues <- NormalizeData(object = my.killi.tissues, normalization.method = "LogNormalize",  scale.factor = 10000)


################################################################################################################################################################
#### 5. Cell cycle regression
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "../../cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

# make into mouse gene names
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

cc.genes.mouse <- firstup(tolower(cc.genes))

# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes   <- cc.genes.mouse[1:43]
g2m.genes <- cc.genes.mouse[44:97]

# Read Blast results for killi 2015 genome to Mouse (Param Singh)
my.nf_mm.homo <- read.table('../../2021-08-19_Killifish_tissue_scRNAseq/Homolos/BestHits_nfur-mmus_1e-3.txt', sep = "\t", header = F)
my.nf_mm.homo$killi_gname <- unlist(lapply(strsplit(my.nf_mm.homo$V1,"|", fixed = TRUE),'[',1))
my.nf_mm.homo$mouse_gname <- unlist(lapply(strsplit(my.nf_mm.homo$V2,"|", fixed = TRUE),'[',1))

# get the cell cycle homologs
s.genes.k  <- my.nf_mm.homo$killi_gname[my.nf_mm.homo$mouse_gname %in% s.genes]
g2m.genes.k <- my.nf_mm.homo$killi_gname[my.nf_mm.homo$mouse_gname %in% g2m.genes]

# Assign Cell-Cycle Scores
my.killi.tissues <- CellCycleScoring(object = my.killi.tissues, s.features = s.genes.k, g2m.features = g2m.genes.k, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x = my.killi.tissues@meta.data)

# write predictions to file
write.table(my.killi.tissues@meta.data, file = paste(Sys.Date(),"Killifish_Tissue_Cohort3_CellCycle_predictions.txt", sep = "_"), sep = "\t", quote = F)

pdf(paste0(Sys.Date(), "_Killifish_Tissue_Cohort3Cohort3_cell_cycle_pie_charts.pdf"), height = 16, width= 10)
par(mfrow=c(3,2))
pie(table(my.killi.tissues@meta.data$Phase[my.killi.tissues@meta.data$Condition == "Female_Kidney_Cohort3"])   , main = "Female_Kidney")
pie(table(my.killi.tissues@meta.data$Phase[my.killi.tissues@meta.data$Condition == "Female_Liver_Cohort3" ])   , main = "Female_Liver" )
pie(table(my.killi.tissues@meta.data$Phase[my.killi.tissues@meta.data$Condition == "Female_Spleen_Cohort3"])   , main = "Female_Spleen")
pie(table(my.killi.tissues@meta.data$Phase[my.killi.tissues@meta.data$Condition == "Male_Kidney_Cohort3"  ])   , main = "Male_Kidney"  )
pie(table(my.killi.tissues@meta.data$Phase[my.killi.tissues@meta.data$Condition == "Male_Liver_Cohort3"   ])   , main = "Male_Liver"   )
pie(table(my.killi.tissues@meta.data$Phase[my.killi.tissues@meta.data$Condition == "Male_Spleen_Cohort3"  ])   , main = "Male_Spleen"  )
par(mfrow=c(1,1))
dev.off()

################################################################################################################################################################
##### 6.Find and remove doublets using doublet finder & Doublet Decon workflow

# https://github.com/chris-mcginnis-ucsf/DoubletFinder

my.killi.tissues <- SCTransform(object = my.killi.tissues, vars.to.regress = c("nFeature_RNA", "percent.mito", "Phase"))
save(my.killi.tissues, file = paste0(Sys.Date(),"_Killi_Tissues_cohort3_Seurat_object_postSCT.RData"))

my.killi.tissues <- RunPCA(my.killi.tissues, npcs = 30)

# Determine the ‘dimensionality’ of the dataset
pdf(paste0(Sys.Date(), "_Killi_Tissue_ElbowPlot.pdf"))
ElbowPlot(my.killi.tissues, ndims = 30)
dev.off()

# run dimensionality reduction
my.killi.tissues <- RunUMAP(my.killi.tissues, dims = 1:30)
my.killi.tissues <- FindNeighbors(my.killi.tissues, dims = 1:30)
my.killi.tissues <- FindClusters(object = my.killi.tissues)

# get cell numbers for each sample, so as to get predicted doublet rate from 10x manual
table(my.killi.tissues$Condition)
# Female_Kidney_Cohort3  Female_Liver_Cohort3 Female_Spleen_Cohort3   Male_Kidney_Cohort3    Male_Liver_Cohort3   Male_Spleen_Cohort3 
# 4649                  2787                  3213                  4405                  2773                  3919 

# Targeted Cell Recovery  # of Cells Loaded	Barcodes Detected	Singlets Multiplets	Multiplet Rate
# 3,000	                         4,950           ~3,000        ~2,900	   ~80	    ~2.4%
# 4,000	                         6,600           ~3,900	       ~3,800	   ~140     ~3.2%
# 5,000	                         8,250           ~4,800        ~4,600	   ~210     ~4.0%
# 6,000                          9,900           ~5,700	       ~5,400	   ~300	    ~4.8%
# 7,000                          11,550	         ~6,600	       ~6,200	   ~400	    ~5.6%
# 8,000                          13,200	         ~7,500	       ~7,000	   ~510	    ~6.4%
# 9,000                          14,850	         ~8,400	       ~7,700	   ~640	    ~7.2%

# assume double the rate for each method, and take the intersection to be more conservative
pred.db.rate <- c(0.040, 0.024, 0.032, 0.040, 0.024, 0.032)

########################################################
#### need to split by 10x sample to make sure to identify real doublets
# will run on one object at a time
cohort3.list <- SplitObject(my.killi.tissues, split.by = "Condition")

########################################################
######## A. Run DoubletFinder
# https://github.com/chris-mcginnis-ucsf/DoubletFinder

# loop over sample
for (i in 1:length(cohort3.list)) {
  
  ## pK Identification (no ground-truth)
  sweep.res.list_killi <- paramSweep_v3(cohort3.list[[i]], PCs = 1:30, sct = TRUE, num.cores	 = 4)
  sweep.stats_killi    <- summarizeSweep(sweep.res.list_killi, GT = FALSE)
  bcmvn_killi          <- find.pK(sweep.stats_killi)
  
  # need some R gymnastics since the Pk is stored as a factor for some reason
  # to get the pK number, need to first convert to character and THEN to numeric
  # numeric first yield row number
  pk.killi <- as.numeric(as.character(bcmvn_killi[as.numeric(bcmvn_killi$pK[bcmvn_killi$BCmetric == max(bcmvn_killi$BCmetric)]),"pK"]))
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(cohort3.list[[i]]@meta.data$seurat_clusters)   ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi       <- round((pred.db.rate[[i]]+ 0.025)                              ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
                          *length(cohort3.list[[i]]@meta.data$orig.ident))        ## Assuming 10x provided doublet formation rate, based on observed cell yield
  
  ## Run DoubletFinder with varying classification stringencies
  cohort3.list[[i]] <- doubletFinder_v3(cohort3.list[[i]], PCs = 1:30, pN = 0.25, pK = pk.killi, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)
  
  # get classification name
  my.DF.res.col <- colnames(cohort3.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort3.list[[i]]@meta.data))]
  
  # rename column to enable subsetting
  colnames(cohort3.list[[i]]@meta.data)[grep("DF.classifications_0.25",colnames(cohort3.list[[i]]@meta.data))] <- "DoubletFinder"
  
}

# run UMAP plots
for (i in 1:length(cohort3.list)) {
  
  pdf(paste(Sys.Date(),"Killifish_Tissue",names(cohort3.list)[i],"Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 5)
  print(DimPlot(cohort3.list[[i]], reduction = "umap", group.by = "DoubletFinder"))
  dev.off()
}

# Remerge the objects post doubletFinder doublet calling
# https://satijalab.org/seurat/articles/merge_vignette.html
killi.singlets.annot.c3 <- merge(cohort3.list[[1]],
                                 y = c(cohort3.list[[2]],
                                       cohort3.list[[3]],
                                       cohort3.list[[4]],
                                       cohort3.list[[5]],
                                       cohort3.list[[6]]),
                                 project = "Killi_Tissues_cohort3")
killi.singlets.annot.c3

# remove pANN columns that are 10xGenomics library lane specific
killi.singlets.annot.c3@meta.data <- killi.singlets.annot.c3@meta.data[,-grep("pANN",colnames(killi.singlets.annot.c3@meta.data))]

########################################################
######## B. Run scds:single cell doublet scoring (hybrid method)
# https://www.bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html
# cxds is based on co-expression of gene pairs and works with absence/presence calls only, 
# bcds uses the full count information and a binary classification approach using artificially generated doublets. 
# cxds_bcds_hybrid combines both approaches

# create scds working object
# convert list to SingleCellExperiment
# https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
cohort3.list.scds        <- lapply(cohort3.list, as.SingleCellExperiment)

# loop over sample
for (i in 1:length(cohort3.list.scds)) {
  
  # Annotate doublets using co-expression based doublet scoring:
  cohort3.list.scds[[i]] <- cxds_bcds_hybrid(cohort3.list.scds[[i]])
  
  # predicted doublet rate
  n.db <- round((pred.db.rate[[i]]+ 0.025)*ncol(cohort3.list.scds[[i]]))                         ## use predicted doublet rate from 10x + 2.5% of possibly non dissociated cells (estimating)
  
  # sort prediction, get top n.db cells
  srt.db.score <- sort(cohort3.list.scds[[i]]$hybrid_score, index.return = T, decreasing = T)
  cohort3.list.scds[[i]]$scds <- "Singlet"
  cohort3.list.scds[[i]]$scds[srt.db.score$ix[1:n.db]] <- "Doublet"
  
}

# run UMAP plots (didn't save)
for (i in 1:length(cohort3.list.scds)) {
  
  p <- plotReducedDim(cohort3.list.scds[[i]], dimred = "UMAP", colour_by = "scds")
  
  pdf(paste(Sys.Date(),"Killifish_Tissue", names(cohort3.list.scds)[i],"scds_UMAP.pdf", sep = "_"), height = 5, width = 5)
  plot(p)
  dev.off()
}

## gate back to doubletFinder annotated Seurat object
killi.singlets.annot.c3@meta.data$scds_hybrid <- NA # initialize

for (i in 1:length(cohort3.list.scds)) {
  
  # for each object compare and move doublet annotations over
  killi.singlets.annot.c3@meta.data[colnames(cohort3.list.scds[[i]]), ]$scds_hybrid <- cohort3.list.scds[[i]]$scds
  
}

########################################################
######## C. Merge and summarize doublet findings

table(killi.singlets.annot.c3@meta.data$DoubletFinder,killi.singlets.annot.c3@meta.data$scds_hybrid)
#          Doublet Singlet
# Doublet     154    1113
# Singlet    1113   19367

killi.singlets.annot.c3@meta.data$DoubletCall <- ifelse( bitOr(killi.singlets.annot.c3@meta.data$DoubletFinder == "Doublet", killi.singlets.annot.c3@meta.data$scds_hybrid == "Doublet") > 0, 
                                                         "Doublet", "Singlet")
table(killi.singlets.annot.c3@meta.data$DoubletCall)
# Doublet Singlet 
#    2380   19367 


# re-run dimensionality reduction for plotting purposes
killi.singlets.annot.c3 <- SCTransform(object = killi.singlets.annot.c3, vars.to.regress = c("nFeature_RNA", "percent.mito", "Phase"))
killi.singlets.annot.c3 <- RunPCA(killi.singlets.annot.c3, npcs = 30)
killi.singlets.annot.c3 <- RunUMAP(killi.singlets.annot.c3, dims = 1:30)

pdf(paste0(Sys.Date(),"_Killifish_Tissues_Cohort3_DoubletCall_UMAP.pdf"))
DimPlot(killi.singlets.annot.c3, reduction = "umap", group.by = "DoubletCall")
dev.off()

# save annotated object
save(killi.singlets.annot.c3, file = paste0(Sys.Date(),"_Killifish_Tissues_cohort3_Seurat_object_with_AnnotatedDoublets.RData"))

### extract/subset only singlets
# save data for singlets df
killi.singlets.c3   <- subset(killi.singlets.annot.c3, subset = DoubletCall %in% "Singlet")  # only keep singlets
killi.singlets.c3
# An object of class Seurat 
# 31560 features across 19367 samples within 2 assays 
# Active assay: SCT (15780 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

table(killi.singlets.c3@meta.data$Condition)
# Female_Kidney_Cohort3  Female_Liver_Cohort3 Female_Spleen_Cohort3   Male_Kidney_Cohort3    Male_Liver_Cohort3   Male_Spleen_Cohort3 
#                  4125                  2523                  2866                  3864                  2512                  3477 

# save filtered/annotated object
save(killi.singlets.c3, file = paste0(Sys.Date(),"_Killifish_Tissues_cohort3_Seurat_object_SINGLETS.RData"))
################################################################################################################################################################


#######################
sink(file = paste(Sys.Date(),"_Killifish_Tissues_Cohort3_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()




