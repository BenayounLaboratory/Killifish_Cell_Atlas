setwd('/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2022-02-01_Killifish_scRNAseq_ATLAS_analyses/Machine_Learning')
options(stringsAsFactors = F)

library(Seurat)
library(caret)
library(randomForest)    # random forest
library(kernlab)         # Radial SVM
library(gbm)             # GBM

library(ggplot2)


source('ML_helper_functions.R') # for balanced accuracy calculation


#####################################################################################################################
#### 1. Load annotated up Seurat Objects
load('../Seurat_merged_analysis/2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData')
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


#####################################################################################################################
#### 2. Prepare data for ML

# for balanced ML training: need to sample equal amounts of cells in each sex/cell type so that no unique cell type drives the model 
# The lowest number per Group in Cohort 1 is for female mast cells (48 cells)
killi.clean.c1@meta.data$BioGroupForSampling <- paste0(killi.clean.c1@meta.data$Annotation_v1, killi.clean.c1@meta.data$Sex)
my.BioGroupForSampling <- unique(killi.clean.c1@meta.data$BioGroupForSampling)

# create object to store info of randomly sampled 40 cells from each group
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
  samp.cells <- sample(cell.names, size = 18)
  
  # extract training and testing
  train.cell.list[[i]] <- samp.cells
  test.cell.list[[i]]  <- setdiff(cell.names,samp.cells)
}

length(unlist(train.cell.list)) # 792

# create subsetted Seurat objects by extracting cells
killi.clean.c1.train   <- subset(killi.clean.c1 , cells = unlist(train.cell.list)) # 38740 features across 792 samples within 3 assays
killi.clean.c1.testing <- subset(killi.clean.c1 , cells = unlist(test.cell.list))  # 38740 features across 40367 samples within 3 assays 

# select most variable BUT robustly expressed genes across cell types
killi.clean.c1.train    <- FindVariableFeatures(killi.clean.c1.train, nfeatures = 5000)
killi.clean.variable    <- VariableFeatures(killi.clean.c1.train)

# get average expression over 14 cell types
Idents(killi.clean.c1.train) <- "Annotation_v1"
killi.clean.c1.train.av      <- AverageExpression(killi.clean.c1.train, assays = "RNA")

# select only genes detected in > 20/22 cell types
my.robust <-  apply(killi.clean.c1.train.av$RNA>0,1,sum) >= 20 # detected in 20 or more cell types
sum(my.robust) # 3043

# genes used as features are noth variable AND robustly expressed
select.genes <- intersect(killi.clean.variable,names(my.robust)[my.robust])
length(select.genes) # 289

# build data matrices for ML (training, testing and validation sets)
train.featmat  <- cbind(killi.clean.c1.train@meta.data[,c("Sex","Tissue", "Annotation_v1" ,"percent.mito")]  ,t(as.matrix(killi.clean.c1.train@assays$SCT@counts[select.genes,])))
test.featmat   <- cbind(killi.clean.c1.testing@meta.data[,c("Sex","Tissue", "Annotation_v1" ,"percent.mito")],t(as.matrix(killi.clean.c1.testing@assays$SCT@counts[select.genes,])))
valid.featmat  <- cbind(killi.clean.c2c3@meta.data[,c("Sex","Tissue", "Annotation_v1" ,"percent.mito")]      ,t(as.matrix(killi.clean.c2c3@assays$SCT@counts[select.genes,])))

dim(train.featmat)
# 792 293
dim(test.featmat)
# 40367   293
dim(valid.featmat)
# 40198   293

# save
save(train.featmat, test.featmat, valid.featmat,
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

# train model with caret train function
my.rf.fit       <- train( Sex ~ .,
                          data       =  train.featmat,
                          method     = "rf",
                          importance = TRUE,
                          trControl  =  my.ctrl.opt.rf,
                          tuneGrid   =  fineGrid.rf,
                          metric     = "balancedAcc")

save(my.rf.fit, file = paste0(Sys.Date(),"_RF_model_Sex.RData"))


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
                            n.minobsinnode    = 50)

# train model with caret train function
my.gbm.fit       <- train( Sex ~ .,
                           data       = train.featmat,
                           method     = "gbm",
                           trControl  = my.ctrl.opt.gbm,
                           tuneGrid   = fineGrid.gbm,
                           metric     = "balancedAcc")

save(my.gbm.fit, file = paste0(Sys.Date(),"_GBM_model_Sex.RData"))



#####################################################################################################################
#### 4. OOB performance analysis
# 
results <- resamples(list("GBM"        = my.gbm.fit    ,
                          "RF"         = my.rf.fit     ))

# summary of model differences
my.model.summaries <- summary(results)

my.10cv.data <- my.model.summaries$values
my.median.acc <- apply(my.10cv.data,2,median)

pdf(paste(Sys.Date(),"Machine_Learning_OOB_Accuracy_boxplots_10CV_beeswarm_RF_GBM.pdf", sep =""), height = 4, width = 6)
par(oma=c(0.1,5,0.1,0.1))
boxplot(my.10cv.data[,order(my.median.acc)], las = 1, ylim = c(0.4,1), horizontal = T, xlab = "10 Fold-CV balanced accuracy", outline = F, col = "purple")
beeswarm::beeswarm(my.10cv.data[,order(my.median.acc)], add = T, horizontal = T, pch = 16)
abline(v = 0.5, col = "red", lty = "dashed")
dev.off()


#####################################################################################################################
#### 5. testing and validation performance analysis

# create output object (for balanced accuracy)
perf.mat           <- matrix(NA,2,4)
colnames(perf.mat) <- c("Training", "Testing", "Validation_c2", "Validation_c3" )
rownames(perf.mat) <- c("RF", "GBM")

# RF
rf.predict.sex.train <- predict(my.rf.fit, train.featmat)
rf.confus.train      <- confusionMatrix(rf.predict.sex.train, as.factor(train.featmat$Sex))
perf.mat["RF", 1]    <- rf.confus.train$byClass[11]

rf.predict.sex.test <- predict(my.rf.fit, test.featmat)
rf.confus.test      <- confusionMatrix(rf.predict.sex.test, as.factor(test.featmat$Sex))
perf.mat["RF", 2]    <- rf.confus.test$byClass[11]

rf.predict.sex.valid <- predict(my.rf.fit, valid.featmat)
rf.confus.valid      <- confusionMatrix(rf.predict.sex.valid, as.factor(valid.featmat$Sex))
# perf.mat["RF", 3]    <- rf.confus.valid$byClass[11]

# get c2 and c3 data
rf.confus.valid.c2      <- confusionMatrix(rf.predict.sex.valid[grep("cohort2",rownames(valid.featmat))], as.factor(valid.featmat$Sex[grep("cohort2",rownames(valid.featmat))]))
rf.confus.valid.c3      <- confusionMatrix(rf.predict.sex.valid[grep("Cohort3",rownames(valid.featmat))], as.factor(valid.featmat$Sex[grep("Cohort3",rownames(valid.featmat))]))
perf.mat["RF", 3]    <- rf.confus.valid.c2$byClass[11]
perf.mat["RF", 4]    <- rf.confus.valid.c3$byClass[11]

# GBM
gbm.predict.sex.train <- predict(my.gbm.fit, train.featmat)
gbm.confus.train      <- confusionMatrix(gbm.predict.sex.train, as.factor(train.featmat$Sex))
perf.mat["GBM", 1]    <- gbm.confus.train$byClass[11]

gbm.predict.sex.test <- predict(my.gbm.fit, test.featmat)
gbm.confus.test      <- confusionMatrix(gbm.predict.sex.test, as.factor(test.featmat$Sex))
perf.mat["GBM", 2]    <- gbm.confus.test$byClass[11]

gbm.predict.sex.valid <- predict(my.gbm.fit, valid.featmat)
gbm.confus.valid      <- confusionMatrix(gbm.predict.sex.valid, as.factor(valid.featmat$Sex))
# perf.mat["GBM", 3]    <- gbm.confus.valid$byClass[11]


# get c2 and c3 data
gbm.confus.valid.c2      <- confusionMatrix(gbm.predict.sex.valid[grep("cohort2",rownames(valid.featmat))], as.factor(valid.featmat$Sex[grep("cohort2",rownames(valid.featmat))]))
gbm.confus.valid.c3      <- confusionMatrix(gbm.predict.sex.valid[grep("Cohort3",rownames(valid.featmat))], as.factor(valid.featmat$Sex[grep("Cohort3",rownames(valid.featmat))]))
perf.mat["GBM", 3]       <- gbm.confus.valid.c2$byClass[11]
perf.mat["GBM", 4]       <- gbm.confus.valid.c3$byClass[11]

pdf(paste(Sys.Date(),"Machine_Learning_Training_Testing_Validation_Accuracy_dotchart_RF_GBM_c2_c3_sep.pdf", sep =""), height = 5, width = 5)
dotchart(perf.mat, xlim = c(0.5,1), pch = 16, pt.cex = 2, main ="ML Sex / Balanced Accuracy")
abline(v = 0.5, col = "red", lty = "dashed")
dev.off()


#####################################################################################################################
#### 6. per cell type performance analysis

cell.types <- sort(unique(killi.combined.filt$Annotation_v1))

############ RF
# parse RF accuracies per cell type
RF.accuracies           <- matrix(NA, length(cell.types), 3)
colnames(RF.accuracies) <- c("RF_Training", "RF_Testing", "RF_Validation")
rownames(RF.accuracies) <- cell.types

# Training
rf.data_and_pred.train <- cbind(train.featmat[,c("Sex","Annotation_v1")],"Preds" = rf.predict.sex.train)

for (i in 1:length(cell.types) ) {
  
  # data in cell type
  data.celltype <- rf.data_and_pred.train[rf.data_and_pred.train$Annotation_v1 == cell.types[i],]
  
  # get balanced accuracy
  RF.accuracies[i,1] <- (confusionMatrix(data.celltype$Preds, as.factor(data.celltype$Sex)))$byClass[11]
  
}

# Testing
rf.data_and_pred.test <- cbind(test.featmat[,c("Sex","Annotation_v1")],"Preds" = rf.predict.sex.test)

for (i in 1:length(cell.types) ) {
  
  # data in cell type
  data.celltype <- rf.data_and_pred.test[rf.data_and_pred.test$Annotation_v1 == cell.types[i],]
  
  # get balanced accuracy
  if (length(unique(data.celltype$Sex)) == 1) {
    RF.accuracies[i,2] <- NA
  } else {
    RF.accuracies[i,2] <- (confusionMatrix(data.celltype$Preds, as.factor(data.celltype$Sex)))$byClass[11]
  }
}


# Validation
rf.data_and_pred.valid <- cbind(valid.featmat[,c("Sex","Annotation_v1")],"Preds" = rf.predict.sex.valid)

for (i in 1:length(cell.types) ) {
  
  # data in cell type
  data.celltype <- rf.data_and_pred.valid[rf.data_and_pred.valid$Annotation_v1 == cell.types[i],]
  
  # get balanced accuracy
  RF.accuracies[i,3] <- (confusionMatrix(data.celltype$Preds, as.factor(data.celltype$Sex)))$byClass[11]
  
}

### Get cell type-wise accuracy heatmap
my.accuracy.colors <- c("white","floralwhite","lightsalmon","indianred1","firebrick3","firebrick","firebrick4")

acc.heat <- rbind(RF.accuracies,
                  c(0.5,0.5,0.5))
rownames(acc.heat)[nrow(acc.heat)] <- c("Random")

pdf(paste0(Sys.Date(),"_Cell_Type_Testing_RF_Balanced_Accuracy_heatmap_CLEAN.pdf"))
pheatmap::pheatmap(acc.heat,
                   cluster_cols = F,
                   cluster_rows = F,
                   col = colorRampPalette(my.accuracy.colors)(50),
                   main = "RF accuracies",
                   cellwidth = 20,
                   cellheight = 20)
dev.off()



############ GBM
# parse GBM accuracies per cell type
GBM.accuracies           <- matrix(NA, length(cell.types), 3)
colnames(GBM.accuracies) <- c("GBM_Training", "GBM_Testing", "GBM_Validation")
rownames(GBM.accuracies) <- cell.types

# Training
gbm.data_and_pred.train <- cbind(train.featmat[,c("Sex","Annotation_v1")],"Preds" = gbm.predict.sex.train)

for (i in 1:length(cell.types) ) {
  
  # data in cell type
  data.celltype <- gbm.data_and_pred.train[gbm.data_and_pred.train$Annotation_v1 == cell.types[i],]
  
  # get balanced accuracy
  GBM.accuracies[i,1] <- (confusionMatrix(data.celltype$Preds, as.factor(data.celltype$Sex)))$byClass[11]
  
}

# Testing
gbm.data_and_pred.test <- cbind(test.featmat[,c("Sex","Annotation_v1")],"Preds" = gbm.predict.sex.test)

for (i in 1:length(cell.types) ) {
  
  # data in cell type
  data.celltype <- gbm.data_and_pred.test[gbm.data_and_pred.test$Annotation_v1 == cell.types[i],]
  
  # get balanced accuracy
  if (length(unique(data.celltype$Sex)) == 1) {
    GBM.accuracies[i,2] <- NA
  } else {
    GBM.accuracies[i,2] <- (confusionMatrix(data.celltype$Preds, as.factor(data.celltype$Sex)))$byClass[11]
  }
  
}


# Validation
gbm.data_and_pred.valid <- cbind(valid.featmat[,c("Sex","Annotation_v1")],"Preds" = gbm.predict.sex.valid)

for (i in 1:length(cell.types) ) {
  
  # data in cell type
  data.celltype <- gbm.data_and_pred.valid[gbm.data_and_pred.valid$Annotation_v1 == cell.types[i],]
  
  # get balanced accuracy
  GBM.accuracies[i,3] <- (confusionMatrix(data.celltype$Preds, as.factor(data.celltype$Sex)))$byClass[11]
  
}

### Get cell type-wise accuracy heatmap
my.accuracy.colors <- c("white","floralwhite","lightsalmon","indianred1","firebrick3","firebrick","firebrick4")

acc.heat <- rbind(GBM.accuracies,
                  c(0.5,0.5,0.5))
rownames(acc.heat)[nrow(acc.heat)] <- c("Random")

pdf(paste0(Sys.Date(),"_Cell_Type_Testing_GBM_Balanced_Accuracy_heatmap_CLEAN.pdf"))
pheatmap::pheatmap(acc.heat,
                   cluster_cols = F,
                   cluster_rows = F,
                   col = colorRampPalette(my.accuracy.colors)(50),
                   main = "GBM accuracies",
                   cellwidth = 20,
                   cellheight = 20)
dev.off()



#####################################################################################################################
#### 7. feature importance
#### https://topepo.github.io/caret/variable-importance.html

############ Caculate/extract feature importance values
# feature importance on tree-based models 
my.rf.Imp    <- varImp(my.rf.fit    , useModel = TRUE , scale = TRUE)$importance
my.gbm.Imp   <- varImp(my.gbm.fit   , useModel = TRUE , scale = TRUE)$importance

my.var.Imps <- data.frame(matrix(0,nrow(my.rf.Imp),5))
colnames(my.var.Imps) <- c("Feature","VarImp_GBM","VarImp_RF","Rank_GBM","Rank_RF")
my.var.Imps$Feature <- rownames(my.rf.Imp)

# parse varImp calculation
for ( i in 1:nrow(my.rf.Imp)) {
  
  my.rf.ix  <- i
  my.gbm.ix <- rownames(my.gbm.Imp) %in% rownames(my.rf.Imp)[i]
  
  my.var.Imps$VarImp_RF[i]    <-  apply(my.rf.Imp,1,mean)[my.rf.ix]
  my.var.Imps$VarImp_GBM[i]   <-  my.gbm.Imp$Overall[my.gbm.ix]
  
}

# now get ranks by decreasing order
my.sort.rf   <- sort(my.var.Imps$VarImp_RF  , index.return = T, decreasing = TRUE)
my.sort.gbm  <- sort(my.var.Imps$VarImp_GBM , index.return = T, decreasing = TRUE)

my.var.Imps$Rank_RF[my.sort.rf$ix]    <- 1:nrow(my.var.Imps)
my.var.Imps$Rank_GBM[my.sort.gbm$ix]  <- 1:nrow(my.var.Imps)

my.var.Imps$RANK_Product               <-  my.var.Imps$Rank_RF * my.var.Imps$Rank_GBM

save(my.var.Imps, file = paste0(Sys.Date(),"_Variable_Importance_Parsing_RF_GBM_CARET_SexId.RData"))
write.table(my.var.Imps,file = paste0(Sys.Date(),"_Variable_Importance_Parsing_RF_GBM_CARET_SexId.txt"), quote = F, row.names = F, sep = "\t")

# feature importance on RF-models
my.rf.varimps.native <- my.rf.fit$finalModel$importance[,3:4]
save(my.rf.varimps.native, file = paste0(Sys.Date(),"_Variable_Importance_RF_native_SexId.RData"))


############ plot top predictors across RF/GBM
my.col.palette.acc <- colorRampPalette(rev(c("#FF9999","indianred1","firebrick1","firebrick3","firebrick4")))(1000)

my.rank.sort <- sort (my.var.Imps$RANK_Product, decreasing = F, index.return = T)
my.var.Imps.sorted <- my.var.Imps[rev(my.rank.sort$ix[1:10]),]

# format for ggplot; initialize data structure separately for RF/GBM, then merge objects
my.rf.results   <- data.frame('Feature_name' =  my.var.Imps.sorted$Feature    ,
                              'VarImp'       =  my.var.Imps.sorted$VarImp_RF ,
                              'Rank'         =  my.var.Imps.sorted$Rank_RF   ,
                              'Model'        = "RF" )

my.gbm.results <- data.frame('Feature_name'  =  my.var.Imps.sorted$Feature   ,
                             'VarImp'       =  my.var.Imps.sorted$VarImp_GBM ,
                             'Rank'         =  my.var.Imps.sorted$Rank_GBM   ,
                             'Model'        = "GBM")

my.complete.results <- rbind(my.rf.results,
                             my.gbm.results)

# will make Importance the size and Rank the color
my.complete.results$Feature_name <- factor(my.complete.results$Feature_name, levels = unique(my.complete.results$Feature_name)) # to preserve the wanted plotting order


pdf(paste(Sys.Date(),"RF_GBM_Models_VarImps_balloon_plot_top_10_SexId.pdf", sep="_"), onefile=F, height = 5, width=6)
my.plot <- ggplot(my.complete.results,aes(x=Model,y=Feature_name,colour=Rank,size=VarImp))+ theme_bw()+ geom_point(shape = 16) 
my.plot <- my.plot + ggtitle("Model Feature Importance") + labs(x = "Model", y = "Feature")
my.plot <- my.plot + scale_colour_gradientn(colours = my.col.palette.acc,space = "Lab", na.value = "grey50", guide = "colourbar")
print(my.plot)
dev.off()  




############ plot violin expression distribution of top predictors
killi.combined.filt <- SetIdent(killi.combined.filt, value = "Annotation_v1") 

# remove special character that create retrieval issues
my.feats <- gsub("`","",my.var.Imps.sorted$Feature)

pdf(paste0(Sys.Date(),"_Seurat_violinPlot_expression_top10_predictors.pdf"), height = 8, width = 12)
VlnPlot(killi.combined.filt, 
        features = rev(my.feats), 
        split.by = "Sex", cols = c("deeppink","deepskyblue"),
        pt.size = 0, assay = 'RNA', stack = T, flip = T)
dev.off()  
# 

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


pdf(paste0(Sys.Date(),"_Seurat_DotPlot_expression_top10_predictors.pdf"), height = 8, width = 10)
DotPlot(killi.combined.filt, features = my.feats, group.by = "Sex", cols = c("lightgrey","firebrick"),   dot.scale = 10) +coord_flip()
dev.off()  

write.table(my.feats, file = paste0(Sys.Date(),"_Names_top10_predictors.txt"), col.names = F, row.names = F, quote = F)


#######################
sink(file = paste(Sys.Date(),"Machine_Learning_R_session_Info_KilliSex.txt", sep =""))
sessionInfo()
sink()



