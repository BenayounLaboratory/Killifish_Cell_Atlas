setwd('/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2022-02-01_Killifish_scRNAseq_ATLAS_analyses/Machine_Learning/ML_per_cell_Type/')
options(stringsAsFactors = F)

library(Seurat)
library(caret)
library(randomForest)    # random forest
library(gbm)             # GBM

library(ggplot2)

library(matrixStats) # for column multiplication


source('../ML_helper_functions.R') # for balanced accuracy calculation

# 2022-02-10
# do ML to see if we can predict biological sex across cell types

# 2023-04-04
# bug fix with Bryan

# 2023-06-15
# change the cell types and methods to respond to reviewer concerns

# 2023-06-19
# Try to train 1 model per cell type

#####################################################################################################################
#### 1. Load annotated up Seurat Objects
load('../../Seurat_merged_analysis/2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData')
killi.combined.filt
# An object of class Seurat 
# 38740 features across 81357 samples within 3 assays 
# Active assay: SCT (18370 features, 5000 variable features)
# 2 other assays present: RNA, integrated
# 2 dimensional reductions calculated: pca, umap

#### For training/testing partition
## We will use cohort 1 vs. cohorts 2 and 3 for training/testing vs. validation set (independent animals)
## This allows for all tissues/sex to be represented in each set of data derived from independent animals/cohorts
killi.clean.c1     <- subset(killi.combined.filt , subset = Batch %in% "Cohort_1")               # 41159
killi.clean.c2c3   <- subset(killi.combined.filt , subset = Batch %in% c("Cohort_2","Cohort_3")) # 40198

# determine cell type/sex distribution
table(killi.clean.c1@meta.data$Sex, killi.clean.c1@meta.data$Annotation_v1)
#    B_Cell_Progenitors B_cells Cholangiocytes Endothelial Erythrocyte_Progenitors Erythrocytes Fibroblasts Hepatocytes Hepatocytes_Efferocytosing
# F                357     829             45         205                    1424         9405          58        1879                       1689
# M                209     482             46         979                     916         9560          93        2711                       2813
# 
#    HSPCs Kidney_distal_tubule Kidney_prox_tubule Lymphoid_progenitors Macrophages Mast_cells Multipotent_progenitors Myeloid_progenitors
# F   267                  103                291                   85         366         50                      60                  27
# M   231                  138                464                   77         550         98                      92                  70
# 
#    Neutrophil_Progenitors Neutrophils NK_T_cells NK_T_progenitor_cells Thrombocytes
# F                    234        1612        452                    18          158
# M                    186        1229        343                    21          237

table(killi.clean.c2c3@meta.data$Sex, killi.clean.c2c3@meta.data$Annotation_v1)
#    B_Cell_Progenitors B_cells Cholangiocytes Endothelial Erythrocyte_Progenitors Erythrocytes Fibroblasts Hepatocytes Hepatocytes_Efferocytosing
# F                654    1524             37         194                    1792         7041          45         757                       2237
# M                188     942             41         204                     712         9445          97        1230                       1975
# 
#    HSPCs Kidney_distal_tubule Kidney_prox_tubule Lymphoid_progenitors Macrophages Mast_cells Multipotent_progenitors Myeloid_progenitors
# F   440                   62                241                  124         528        378                     108                  95
# M   226                    2                 24                   45         932        701                      80                 122
# 
#    Neutrophil_Progenitors Neutrophils NK_T_cells NK_T_progenitor_cells Thrombocytes
# F                    435        2687        467                    36          198
# M                    218        1849        641                    80          364


############### ONLY keep cell types with at least 100 cells per sex in cohort 1
# nb. of cells per cluster-sample
cell.per.samp.tab <- t(table(killi.clean.c1@meta.data$Sex, killi.clean.c1@meta.data$Annotation_v1))

# cell types with at least 100 cells from every each sex/cohort sample
celltype.qc <- rownames(cell.per.samp.tab)[apply(cell.per.samp.tab > 100, 1, sum) == 2]

# filter
killi.clean.c1     <- subset(killi.clean.c1   , subset = Annotation_v1 %in% celltype.qc)   # 40319
killi.clean.c2c3   <- subset(killi.clean.c2c3 , subset = Annotation_v1 %in% celltype.qc)   # 38209


#####################################################################################################################
#### 2. Prepare data for ML

# for ML training: will train a RF model for each cell type, need to sample cells in each sex/cell type
# get 100 cells per group
killi.clean.c1@meta.data$BioGroupForSampling <- paste0(killi.clean.c1@meta.data$Annotation_v1, killi.clean.c1@meta.data$Sex)
my.BioGroupForSampling <- unique(killi.clean.c1@meta.data$BioGroupForSampling)

# create object to store info of randomly sampled cells from each group
train.cell.list        <- vector(mode = "list", length = length(my.BioGroupForSampling))
names(train.cell.list) <- my.BioGroupForSampling

# create object to store info of the remaining cells from each group
test.cell.list        <- vector(mode = "list", length = length(my.BioGroupForSampling))
names(test.cell.list) <- my.BioGroupForSampling

# set seed for reproducibility
set.seed(123456789)

for (i in 1:length(my.BioGroupForSampling)){
  # identify cells from the group
  cell.names <- rownames(killi.clean.c1@meta.data)[killi.clean.c1@meta.data$BioGroupForSampling %in% my.BioGroupForSampling[i]]
  
  # sample 18 cells from group[smallest group in cohort 1]
  samp.cells <- sample(cell.names, size = 100)
  
  # extract training and testing
  train.cell.list[[i]] <- samp.cells
  test.cell.list[[i]]  <- setdiff(cell.names,samp.cells)
}

length(unlist(train.cell.list)) # 3000

# create subsetted Seurat objects by extracting cells
killi.clean.c1.train   <- subset(killi.clean.c1 , cells = unlist(train.cell.list)) # 38740 features across 3000 samples within 3 assays
killi.clean.c1.testing <- subset(killi.clean.c1 , cells = unlist(test.cell.list))  # 38740 features across 37319 samples within 3 assays

# select most variable BUT robustly expressed genes across cell types
killi.clean.c1.train    <- FindVariableFeatures(killi.clean.c1.train, nfeatures = 5000)
killi.clean.variable    <- VariableFeatures(killi.clean.c1.train)

# get average expression over 15 cell types
Idents(killi.clean.c1.train) <- "Annotation_v1"
killi.clean.c1.train.av      <- AverageExpression(killi.clean.c1.train, assays = "RNA")

# select only genes detected in all cell types
my.robust <-  apply(killi.clean.c1.train.av$RNA>0,1,sum) == 15 # detected in all 15 cell types
sum(my.robust) # 5417

# genes used as features are both variable AND robustly expressed
select.genes <- intersect(killi.clean.variable,names(my.robust)[my.robust])
length(select.genes) # 495

# build data matrices for ML (training, testing and validation sets)
train.featmat  <- cbind(killi.clean.c1.train@meta.data[,c("Sex","Tissue", "Annotation_v1" ,"percent.mito")]  ,t(as.matrix(killi.clean.c1.train@assays$SCT@counts[select.genes,])))
test.featmat   <- cbind(killi.clean.c1.testing@meta.data[,c("Sex","Tissue", "Annotation_v1" ,"percent.mito")],t(as.matrix(killi.clean.c1.testing@assays$SCT@counts[select.genes,])))
valid.featmat  <- cbind(killi.clean.c2c3@meta.data[,c("Sex","Tissue", "Annotation_v1" ,"percent.mito")]      ,t(as.matrix(killi.clean.c2c3@assays$SCT@counts[select.genes,])))

dim(train.featmat)
# 3000  499
dim(test.featmat)
# 37319   499
dim(valid.featmat)
# 38209   499

# build lists of matrices for each cell type for ML (training, testing and validation sets)
train.list             <- vector(mode = "list", length = length(celltype.qc))
names(train.list)      <- celltype.qc

test.list             <- vector(mode = "list", length = length(celltype.qc))
names(test.list)      <- celltype.qc

valid.list             <- vector(mode = "list", length = length(celltype.qc))
names(valid.list)      <- celltype.qc

for (i in 1:length(celltype.qc)) {
  
  train.list[[i]] <- train.featmat[train.featmat$Annotation_v1 %in% celltype.qc[i],]
  test.list[[i]]  <- test.featmat [test.featmat$Annotation_v1  %in% celltype.qc[i],]
  valid.list[[i]] <- valid.featmat[valid.featmat$Annotation_v1 %in% celltype.qc[i],]
  
}

# save
save(train.list, test.list, valid.list,
     file = paste0(Sys.Date(),"_Training_Testing_Validation_Data_Killi_Atlas.RData"))

#####################################################################################################################
#### 3. Run  ML

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%% A. run RF model

# set seed for reproducibility
set.seed(123456789)

# use 10-fold cross-validation to build the model
my.ctrl.opt.rf <- trainControl(method          = "cv",
                               number          = 10,
                               allowParallel   = TRUE,
                               verbose         = F,
                               summaryFunction = getTwoClassBalancedAccuracy,
                               classProbs      = TRUE)

fineGrid.rf <- expand.grid(mtry = seq(10,250,20))

# create vector to receive trained models
my.rf.fit.list        <- vector(mode = "list", length = length(celltype.qc))
names(my.rf.fit.list) <- celltype.qc

# Loop over cell types
for (i in 1:length(celltype.qc)) {
  # train model with caret train function
  my.rf.fit.list[[i]]       <- train( Sex ~ .,
                                      data       =  train.list[[i]][,-c(2,3)], # -3 to remove cell type since homogeneous, 2 tissue since creates issue
                                      method     = "rf",
                                      importance = TRUE,
                                      trControl  =  my.ctrl.opt.rf,
                                      tuneGrid   =  fineGrid.rf,
                                      metric     = "balancedAcc")
}

save(my.rf.fit.list, file = paste0(Sys.Date(),"_RF_models_per_CellType_Sex.RData"))


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%% B. run GBM model
# set seed for reproducibility
set.seed(123456789)

# use 10-fold cross-validation to build the model
my.ctrl.opt.gbm <- trainControl(method          = "cv",
                                number          = 10,
                                allowParallel   = TRUE,
                                verbose         = F,
                                summaryFunction = getTwoClassBalancedAccuracy,
                                classProbs      = TRUE)

fineGrid.gbm <- expand.grid(n.trees           = seq(2500,10000,2500),
                            interaction.depth = seq(6,18,3),
                            shrinkage         = 0.001,
                            n.minobsinnode    = 10)

# create vector to receive trained models
my.gbm.fit.list        <- vector(mode = "list", length = length(celltype.qc))
names(my.gbm.fit.list) <- celltype.qc

# Loop over cell types
for (i in 1:length(celltype.qc)) {
  # train model with caret train function
  my.gbm.fit.list[[i]]       <- train( Sex ~ .,
                                       data       = train.list[[i]][,-c(2,3)], # -3 to remove cell type since homogeneous, 2 tissue since creates issue
                                       method     = "gbm",
                                       trControl  = my.ctrl.opt.gbm,
                                       tuneGrid   = fineGrid.gbm,
                                       metric     = "balancedAcc")
}

save(my.gbm.fit.list, file = paste0(Sys.Date(),"_GBM_models_per_CellType_Sex.RData"))



#####################################################################################################################
#### 4. OOB performance analysis

# Grab 10-fold CV OOB training accuracies

####### RF
results.rf  <- resamples(my.rf.fit.list)

# summary of model differences
my.model.summaries.rf <- summary(results.rf)

my.10cv.data.rf <- my.model.summaries.rf$values
# my.median.acc.rf <- apply(my.10cv.data.rf,2,median)

pdf(paste(Sys.Date(),"Machine_Learning_OOB_Accuracy_boxplots_10CV_beeswarm_RF_over_cellTypes.pdf", sep =""), height = 7, width = 8)
par(oma=c(0.1,9,0.1,0.1))
boxplot(rev(my.10cv.data.rf), las = 1, ylim = c(0.4,1), horizontal = T, xlab = "10 Fold-CV balanced accuracy", outline = F, col = "purple", main = "RF")
beeswarm::beeswarm(rev(my.10cv.data.rf), add = T, horizontal = T, pch = 16)
abline(v = 0.5, col = "red", lty = "dashed")
dev.off()


####### GBM
results.gbm  <- resamples(my.gbm.fit.list)

# summary of model differences
my.model.summaries.gbm <- summary(results.gbm)

my.10cv.data.gbm <- my.model.summaries.gbm$values
# my.median.acc.gbm <- apply(my.10cv.data.gbm,2,median)

pdf(paste(Sys.Date(),"Machine_Learning_OOB_Accuracy_boxplots_10CV_beeswarm_GBM_over_cellTypes.pdf", sep =""), height = 7, width = 8)
par(oma=c(0.1,9,0.1,0.1))
boxplot(rev(my.10cv.data.gbm), las = 1, ylim = c(0.4,1), horizontal = T, xlab = "10 Fold-CV balanced accuracy", outline = F, col = "purple", main = "GBM")
beeswarm::beeswarm(rev(my.10cv.data.gbm), add = T, horizontal = T, pch = 16)
abline(v = 0.5, col = "red", lty = "dashed")
dev.off()


#####################################################################################################################
#### 5. testing and validation performance analysis

######################################################################
# RF
# create output object (for balanced accuracy)
perf.mat.rf           <- matrix(NA,length(celltype.qc),3)
rownames(perf.mat.rf) <- celltype.qc
colnames(perf.mat.rf) <- c("RF_Testing", "RF_Validation_c2", "RF_Validation_c3")

# loop over cell types for models and data
for (i in 1:length(celltype.qc)) {
  
  # use cell type appropriate model
  rf.predict.sex.test                         <- predict(my.rf.fit.list[[i]], test.list[[i]])
  rf.confus.test                              <- confusionMatrix(rf.predict.sex.test, as.factor(test.list[[i]]$Sex))
  perf.mat.rf[celltype.qc[i], "RF_Testing"]   <- rf.confus.test$byClass[11]
  
  # get validation data
  rf.predict.sex.valid                              <- predict(my.rf.fit.list[[i]], (valid.list[[i]]))

  # since c2 is missing some tissues/sex combo, it will throw error if only one sex present
  # in that case, get an NA result
  if(length(unique(valid.list[[i]]$Sex[grep("cohort2",rownames(valid.list[[i]]))])) > 1) {
    rf.confus.valid.c2                                <- confusionMatrix(rf.predict.sex.valid[grep("cohort2",rownames(valid.list[[i]]))], as.factor(valid.list[[i]]$Sex[grep("cohort2",rownames(valid.list[[i]]))]))
    perf.mat.rf[celltype.qc[i], "RF_Validation_c2"]   <- rf.confus.valid.c2$byClass[11]
  } else {
    perf.mat.rf[celltype.qc[i], "RF_Validation_c2"]   <- NA
  }

  rf.confus.valid.c3                                <- confusionMatrix(rf.predict.sex.valid[grep("Cohort3",rownames(valid.list[[i]]))], as.factor(valid.list[[i]]$Sex[grep("Cohort3",rownames(valid.list[[i]]))]))
  perf.mat.rf[celltype.qc[i], "RF_Validation_c3"]   <- rf.confus.valid.c3$byClass[11]
  
}

pdf(paste(Sys.Date(),"Machine_Learning_Testing_Validation_Accuracy_barplot_RF.pdf", sep =""), height = 6, width = 6.5)
par(oma=c(0.1,7,0.1,0.1))
barplot(t(perf.mat.rf[rev(celltype.qc),3:1]), beside = T, horiz = T, las = 1, 
        col = c("steelblue2", "royalblue4", "gray48"), xlim = c(0.4,1), xpd = F,
        xlab = "ML Sex / Balanced Accuracy", main = "RF")
box()
abline(v = 0.5, col = "red", lty = "dashed")

legend("topright", c("Testing", "Validation (Cohort 2)", "Validation (Cohort 3)"), fill = c("gray48","royalblue4", "steelblue2" ), bty = 'n', cex = 0.7)
dev.off()


######################################################################
# GBM

# create output object (for balanced accuracy)
perf.mat.gbm           <- matrix(NA,length(celltype.qc),3)
rownames(perf.mat.gbm) <- celltype.qc
colnames(perf.mat.gbm) <- c("GBM_Testing", "GBM_Validation_c2", "GBM_Validation_c3")

# loop over cell types for models and data
for (i in 1:length(celltype.qc)) {
  
  # use cell type appropriate model
  gbm.predict.sex.test                         <- predict(my.gbm.fit.list[[i]], test.list[[i]])
  gbm.confus.test                              <- confusionMatrix(gbm.predict.sex.test, as.factor(test.list[[i]]$Sex))
  perf.mat.gbm[celltype.qc[i], "GBM_Testing"]   <- gbm.confus.test$byClass[11]
  
  # get validation data
  gbm.predict.sex.valid                              <- predict(my.gbm.fit.list[[i]], (valid.list[[i]]))
  
  # since c2 is missing some tissues/sex combo, it will throw error if only one sex present
  # in that case, get an NA result
  if(length(unique(valid.list[[i]]$Sex[grep("cohort2",rownames(valid.list[[i]]))])) > 1) {
    gbm.confus.valid.c2                                <- confusionMatrix(gbm.predict.sex.valid[grep("cohort2",rownames(valid.list[[i]]))], as.factor(valid.list[[i]]$Sex[grep("cohort2",rownames(valid.list[[i]]))]))
    perf.mat.gbm[celltype.qc[i], "GBM_Validation_c2"]   <- gbm.confus.valid.c2$byClass[11]
  } else {
    perf.mat.gbm[celltype.qc[i], "GBM_Validation_c2"]   <- NA
  }
  
  gbm.confus.valid.c3                                <- confusionMatrix(gbm.predict.sex.valid[grep("Cohort3",rownames(valid.list[[i]]))], as.factor(valid.list[[i]]$Sex[grep("Cohort3",rownames(valid.list[[i]]))]))
  perf.mat.gbm[celltype.qc[i], "GBM_Validation_c3"]   <- gbm.confus.valid.c3$byClass[11]
  
}

pdf(paste(Sys.Date(),"Machine_Learning_Testing_Validation_Accuracy_barplot_GBM.pdf", sep =""), height = 6, width = 6.5)
par(oma=c(0.1,7,0.1,0.1))
barplot(t(perf.mat.gbm[rev(celltype.qc),3:1]), beside = T, horiz = T, las = 1, 
        col = c("steelblue2", "royalblue4", "gray48"), xlim = c(0.4,1), xpd = F,
        xlab = "ML Sex / Balanced Accuracy", main = "GBM")
box()
abline(v = 0.5, col = "red", lty = "dashed")

legend("topright", c("Testing", "Validation (Cohort 2)", "Validation (Cohort 3)"), fill = c("gray48","royalblue4", "steelblue2" ), bty = 'n', cex = 0.7)
dev.off()


#####################################################################################################################
#### 6. feature importance
#### https://topepo.github.io/caret/variable-importance.html

# feats for improtance (-3 for the Sex, and tissue and cell type variables)
nfeats.for.imp <- ncol(test.list[[1]]) -3

############ Caculate/extract feature importance values
# Create df to receive var imps
#### RF
my.var.Imps.RF             <- data.frame(matrix(0,nfeats.for.imp,length((my.rf.fit.list))))
colnames(my.var.Imps.RF)   <- names(my.rf.fit.list)
rownames(my.var.Imps.RF)   <- colnames(train.list[[i]])[-c(1:3)]

for (i in 1:length(my.rf.fit.list)){
  # feature importance on tree-based models
  my.var.Imps.RF[,i] <- varImp(my.rf.fit.list[[i]]    , useModel = TRUE , scale = TRUE)$importance[,1]
}

#### GBM
my.var.Imps.GBM             <- data.frame(matrix(0,nfeats.for.imp,length((my.gbm.fit.list))))
colnames(my.var.Imps.GBM)   <- names(my.gbm.fit.list)
rownames(my.var.Imps.GBM)   <- colnames(train.list[[i]])[-c(1:3)]

for (i in 1:length(my.gbm.fit.list)){
  # feature importance on tree-based models
  my.var.Imps.GBM[,i] <- varImp(my.gbm.fit.list[[i]]    , useModel = TRUE , scale = TRUE)$importance[,1]
}



############ Caculate ranks of VarImp
# Since high importance is better, and caret maxes out at 100, do substraction get low rank to be better
Rank.varImps.RF  <- apply(100-my.var.Imps.RF,2,rank)
Rank.varImps.GBM <- apply(100-my.var.Imps.GBM,2,rank)

my.var.Imps.RF$RANK_Product  <- rowProds(as.matrix(Rank.varImps.RF))
my.var.Imps.GBM$RANK_Product <- rowProds(as.matrix(Rank.varImps.GBM))

# Save feature importance and rank to files
write.table(my.var.Imps.RF, file = paste0(Sys.Date(),"_Variable_Importance_Parsing_RF_perCellType_CARET.txt") , quote = F, row.names = T, sep = "\t")
write.table(my.var.Imps.GBM,file = paste0(Sys.Date(),"_Variable_Importance_Parsing_GBM_perCellType_CARET.txt"), quote = F, row.names = T, sep = "\t")

# for comb
my.comb.VarImps <- merge(my.var.Imps.RF, my.var.Imps.GBM, by = "row.names", suffixes = c("_RF","_GBM"))
colnames(my.comb.VarImps)[1] <- "Feature"

# Get combined Rank product
my.comb.VarImps$Combined_RANK_Product <- my.comb.VarImps$RANK_Product_RF * my.comb.VarImps$RANK_Product_GBM

# Save combined feature importance and rank to files
write.table(my.comb.VarImps, file = paste0(Sys.Date(),"_Variable_Importance_Parsing_perCellType_CARET_RF_GBM_summary.txt") , quote = F, row.names = F, sep = "\t")

save(my.comb.VarImps, file = paste0(Sys.Date(),"_CombinedVarImps.RData"))


############ plot top predictors across RF/GBM
my.col.palette.acc <- colorRampPalette(rev(c("#FF9999","indianred1","firebrick1","firebrick3","firebrick4")))(1000)

my.rank.sort       <- sort (my.comb.VarImps$Combined_RANK_Product, decreasing = F, index.return = T)
my.var.Imps.sorted <- my.comb.VarImps[rev(my.rank.sort$ix[1:10]),]


# format for ggplot; initialize data structure separately for RF/GBM, then merge objects
my.rf.results   <- data.frame('Feature_name'        =  my.var.Imps.sorted$Feature    ,
                              'Median_VarImp'       =  apply(my.var.Imps.sorted[,2:16],1,median) ,
                              'RankProd'            =  my.var.Imps.sorted$RANK_Product_RF   ,
                              'Algorithm'           = "RF" )

my.gbm.results <- data.frame('Feature_name'        =  my.var.Imps.sorted$Feature   ,
                             'Median_VarImp'       =  apply(my.var.Imps.sorted[,18:32],1,median),
                             'RankProd'            =  my.var.Imps.sorted$RANK_Product_GBM   ,
                             'Algorithm'           = "GBM")

my.complete.results <- rbind(my.rf.results,
                             my.gbm.results)

# will make Importance the size and Rank the color
my.complete.results$Feature_name <- factor(my.complete.results$Feature_name, levels = unique(my.complete.results$Feature_name)) # to preserve the wanted plotting order


pdf(paste(Sys.Date(),"RF_GBM_Models_VarImps_balloon_plot_top_10_SexId_CellTypeMedian.pdf", sep="_"), onefile=F, height = 5, width=6)
my.plot <- ggplot(my.complete.results,aes(x=Algorithm,y=Feature_name,colour=RankProd,size=Median_VarImp))+ theme_bw()+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("Model Feature Importance") + labs(x = "Algorithm", y = "Feature")
my.plot <- my.plot + scale_colour_gradientn(colours = my.col.palette.acc,space = "Lab", na.value = "grey50", guide = "colourbar")
print(my.plot)
dev.off()


############ plot violin expression distribution of top predictors
load("2023-06-22_CombinedVarImps.RData")
load('../../Seurat_merged_analysis/2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData')

killi.combined.filt <- SetIdent(killi.combined.filt, value = "Annotation_v1")
DefaultAssay(killi.combined.filt) <- "RNA"

### top gene wise
pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_LOC107373896_ncfem1.pdf"), height = 4, width = 12)
VlnPlot(killi.combined.filt,
        features = "LOC107373896",
        split.by = "Sex", cols = c("deeppink","deepskyblue"),
        pt.size = 0, assay = 'RNA', flip = T)
dev.off()

pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_LOC107372823_hsp70.pdf"), height = 4, width = 12)
VlnPlot(killi.combined.filt,
        features = "LOC107372823",
        split.by = "Sex", cols = c("deeppink","deepskyblue"),
        pt.size = 0, assay = 'RNA', flip = T)
dev.off()

pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_LINE−R2.pdf"), height = 4, width = 12)
VlnPlot(killi.combined.filt,
        features = "NotFur1-GapFilledScaffold-375-227342-LINE-R2",
        split.by = "Sex", cols = c("deeppink","deepskyblue"),
        pt.size = 0, assay = 'RNA', flip = T)
dev.off()


#############################
# remove special character that create retrieval issues
my.feats <- gsub("`","",my.var.Imps.sorted$Feature)

# code throws out conflict due to SCT/integrated, remove for plotting
killi.combined.filt[['SCT']] <- NULL
killi.combined.filt[['integrated']] <- NULL
DefaultAssay(killi.combined.filt) <- "RNA"

pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_top10_predictors_CellTypeModel.pdf"), height = 8, width = 12)
VlnPlot(killi.combined.filt,
        features = rev(my.feats),
        split.by = "Sex", cols = c("deeppink","deepskyblue"),
        pt.size = 0, assay = 'RNA', stack = T, flip = T)
dev.off()


# get average expression over Sex
Idents(killi.combined.filt) <- "Sex"
killi.combined.filt.av      <- AverageExpression(killi.combined.filt, assays = "RNA")

sort.feats <- sort(killi.combined.filt.av$RNA[my.feats,]$F - killi.combined.filt.av$RNA[my.feats,]$M, index.return = T)


pdf(paste0(Sys.Date(),"_Seurat_DotPlot_expression_top10_predictors_CellTypeModel.pdf"), height = 8, width = 10)
DotPlot(killi.combined.filt, features = paste0("rna_",my.feats), group.by = "Sex", cols = c("lightgrey","firebrick"),   dot.scale = 10, assay = 'RNA') + coord_flip()
dev.off()

pdf(paste0(Sys.Date(),"_Seurat_DotPlot_expression_top10_predictors_CellTypeModel_SORTED.pdf"), height = 8, width = 10)
DotPlot(killi.combined.filt, features = paste0("rna_",my.feats)[sort.feats$ix], group.by = "Sex", cols = c("lightgrey","firebrick"),   dot.scale = 10, assay = 'RNA') + coord_flip()
dev.off()

write.table(my.feats, file = paste0(Sys.Date(),"_Names_top10_predictors.txt"), col.names = F, row.names = F, quote = F)


#######################
sink(file = paste(Sys.Date(),"Machine_Learning_R_session_Info_KilliSex_Per_Cell_Type.txt", sep =""))
sessionInfo()
sink()



