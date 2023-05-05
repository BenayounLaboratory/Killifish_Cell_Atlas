# Script was adapted from CSCDRNA Example Scripts
# https://github.com/empiricalbayes/CSCDRNA

# Setting the working directory
setwd("/scratch2/alanxu/liver_deconv")

# Loading necessary libraries
library(Seurat)
library(CSCDRNA)
library(Biobase)

load("liver_only_for_deconvolution.Rdata")

#Build ExpressionSet with bulk data.
bulk.df <- read.table("liver.txt", sep = "\t", header = TRUE)
bulk.df <- bulk.df[,-(2:6)]
bulk.df$Geneid <- gsub("gene-","",bulk.df$Geneid,fixed = TRUE)
rownames(bulk.df) <- bulk.df$Geneid
bulk.df <- bulk.df[,-1]
bulk.mat <- data.matrix(bulk.df)

bulk.eset <- Biobase::ExpressionSet(assayData = bulk.mat)

#Build ExpressionSet with single-cell data.
sc.counts.matrix=data.matrix(onlyliver@assays$RNA@counts)

individual.labels=onlyliver@meta.data$Condition

cell.type.labels=onlyliver@meta.data$Annotation_v1

sample.ids <- colnames(sc.counts.matrix)

#individual.labels and cell.types should be in the same order as in sample.ids.

sc.pheno <- data.frame(check.names=FALSE, check.rows=FALSE,
                       stringsAsFactors=FALSE,row.names=sample.ids,
                       SubjectName=individual.labels,cellType=cell.type.labels)
sc.meta <- data.frame(labelDescription=c("SubjectName","cellType"),
                      row.names=c("SubjectName","cellType"))
sc.pdata <- new("AnnotatedDataFrame",data=sc.pheno, varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData=sc.counts.matrix,phenoData=sc.pdata)

analysis <- CSCD(bulk.eset=bulk.eset,sc.eset= sc.eset,
                 min.p=0.3,markers=NULL,cell.types="cellType",
                 subj.names="SubjectName",verbose=TRUE)

write.table(analysis$bulk.props,file = "liver_deconv.txt", quote = FALSE)
save(analysis, file = "analysis.Rdata")