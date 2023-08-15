#Quantifying Efferocytosis, ncFem1 expression, and ncFem1 localization

#setwd
setwd("/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/CODE/12_RNAScope_Quantification")

#load libraries
library(readxl)
library(reshape2)
library(beeswarm)
library(outliers)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

#1. Quantifying hepatocyte efferocytosis

#import counts
table <- read.csv("input/Efferocytosis_cell_type_quant.xlsx")
table <- table[,c(1,5:7)]
colnames(table) <- c("Sex","Percentage_Erythrocytes", "Percentage_Hepatocytes", "Percentage_Efferocytosing_Hepatocytes")

table$Sex <- gsub("_slice.+", "", table$Sex)

#take the median values for each group
table_diff <- table %>% group_by(Sex) %>% 
  summarise(med_eryth=median(Percentage_Erythrocytes),
            med_hep= median(Percentage_Hepatocytes),
            med_effhep= median(Percentage_Efferocytosing_Hepatocytes),
            .groups = 'drop') %>%
  as.data.frame()

#change colnames
colnames(table_diff)[c(2:4)] <- c("Erythrocytes", "Hepatocytes", "Hepatocytes_Efferocytosing")

#rearrange scores for plotting
aggregate <- melt(table_diff)
aggregate[is.na(aggregate)] <- 0
aggregate$Sex <- gsub("\\d", "", aggregate$Sex)
aggregate$Sex <- gsub("eff_", "", aggregate$Sex)

#export table
write.table(aggregate, file = paste0(Sys.Date(),  "Cell_type_proportions_RNAScope_livers.txt"), quote = F, sep = "\t", row.names = F)

#full plot
pdf(paste(Sys.Date(), "eff_quant.pdf", sep = "_"), height = 5, width= 5)
ggplot(aggregate, aes(x=variable, y=value)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  geom_jitter(aes(color=Sex), position = position_jitter(height = .1, width = .1), cex = 1) +
  scale_color_manual(values = c('female' = 'Deeppink', 'male' = 'Deepskyblue'))
dev.off()
###

#only efferocytosing hepatocytes plot
eff_only <- subset(aggregate, aggregate$variable == "Hepatocytes_Efferocytosing")
colnames(eff_only) <- c("Sex", "Cell_Type", "Percentage")

#plot
pdf(paste(Sys.Date(), "eff_quant_only.pdf", sep = "_"), height = 5, width= 5)
ggplot(eff_only, aes(x=Cell_Type, y=Percentage)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  geom_jitter(aes(color=Sex), position = position_jitter(height = .1, width = .1), cex = 1) +
  scale_color_manual(values = c('female' = 'Deeppink', 'male' = 'Deepskyblue'))
dev.off()

#plot with both sexes included
pdf(paste(Sys.Date(), "eff_quant_only_sex_sep.pdf", sep = "_"), height = 5, width= 5)
ggplot(eff_only, aes(x=Cell_Type, y=Percentage, fill = Sex)) +
  scale_fill_manual(values = c("deeppink",  "deepskyblue")) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  geom_point(position = position_dodge(width = .75)) +
  stat_compare_means(method = "t.test" )
dev.off()


#2. Quantifying differential ncFem1 expression by sex
# read excel sheet
my.data <- read_xlsx("input/ncFem1_quantification.xlsx", sheet = 1)

# Do stats
dapi.pval  <- wilcox.test(my.data$ncFem1_norm_to_M_DAPI[my.data$Sex   == "F"],
                          my.data$ncFem1_norm_to_M_DAPI[my.data$Sex == "M"])
ppib.pval  <- wilcox.test(my.data$ncFem1_norm_to_M_ppIB[my.data$Sex   == "F"],
                          my.data$ncFem1_norm_to_M_ppIB[my.data$Sex == "M"])

my.data.plotting <- list("dapi_normalized_F" = my.data$ncFem1_norm_to_M_DAPI[my.data$Sex   == "F"],
                         "dapi_normalized_M" = my.data$ncFem1_norm_to_M_DAPI[my.data$Sex   == "M"],
                         "ppib_normalized_F" = my.data$ncFem1_norm_to_M_ppIB[my.data$Sex   == "F"],
                         "ppib_normalized_M" = my.data$ncFem1_norm_to_M_ppIB[my.data$Sex   == "M"])


# make data.frame for plotting
pdf(paste0(Sys.Date(),"_both_housekeepers_normalized_ncfem1_RNAScope_signal.pdf"), height = 5, width = 6)
boxplot( rev(my.data.plotting), outline = F, las = 1,
         col = c("deepskyblue","deeppink"), 
         ylab = 'Ratio of ncFem1 area to housekeeper area', main = "ncFem1:housekeeper", ylim = c(0,5),
         xlab = "", horizontal=F)
beeswarm(rev(my.data.plotting), add = T, pch = 16, cex = 1.2, horizontal=F)
text(1.5, 4.5, signif(dapi.pval$p.value,2))
text(3.5, 4.5,signif(ppib.pval$p.value,2))
abline(h = 1, col = "red", lty = "dashed")
dev.off()


#use only ppib as a housekeeper
my.data.plotting.ppib <- list("ppib_normalized_F" = my.data$ncFem1_norm_to_M_ppIB[my.data$Sex   == "F"],
                         "ppib_normalized_M" = my.data$ncFem1_norm_to_M_ppIB[my.data$Sex   == "M"])


pdf(paste0(Sys.Date(),"_ppib_normalized_ncfem1_RNAScope_signal.pdf"), height = 5, width = 2)
boxplot( rev(my.data.plotting.ppib), outline = F, las = 1,
         col = c("deepskyblue","deeppink"), 
         ylab = 'Ratio of ncFem1 area to housekeeper area', main = "ncFem1:housekeeper", ylim = c(0,5),
         xlab = "", horizontal=F)
beeswarm(rev(my.data.plotting.ppib), add = T, pch = 16, cex = 1.2, horizontal=F)
text(1.5, 4.5,signif(ppib.pval$p.value,2))
abline(h = 1, col = "red", lty = "dashed")
dev.off()


#3. Nuclear ncFem1 localization

#import data
ncfem <- read.csv("input/ncFem1_localization.csv")

ncfem$Image <- gsub(".+-", "", ncfem$Image)

#keep mean values for each channel for nuclei and cytoplasm; reports for each cell
relevant <- ncfem[,c("Image", "Nucleus..Channel.2.mean", "Nucleus..Channel.3.mean", "Cytoplasm..Channel.2.mean", "Cytoplasm..Channel.3.mean")]

#take nuc:cytoplasm ratio
relevant$ppib_ratio <- relevant$Nucleus..Channel.3.mean/relevant$Cytoplasm..Channel.3.mean
relevant$nc_ratio <- relevant$Nucleus..Channel.2.mean/relevant$Cytoplasm..Channel.2.mean

#now take only those values
relevant <- relevant[,c(1,6,7)]

#take mean values for each image
location <- relevant %>% group_by(Image) %>% 
  summarise(mean_nc=mean(nc_ratio),
            mean_ppib= mean(ppib_ratio),
            .groups = 'drop') %>%
  as.data.frame()

#take median value per group
location$biorep <- location$Image
location$biorep <- gsub("_slice.+","", location$biorep)

location_biorep <- location %>% group_by(biorep) %>% 
  summarise(med_ppib=median(mean_ppib),
            med_nc= median(mean_nc),
            .groups = 'drop') %>%
  as.data.frame()

location_biorep$sex <- c(rep("Female", 5), rep("Male", 5))

#rearrange scores for plotting
aggregate <- melt(location_biorep)

pdf(paste(Sys.Date(), "ncfem1_to_ppib_nuclear_signal.pdf", sep = "_"), height = 5, width = 5)
ggplot(aggregate, aes(x=variable, y=value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  scale_fill_manual(values = c("dark red",  "yellow")) +
  geom_jitter(aes(color=sex), position = position_jitter(height = .1, width = .1), cex = 1) +
  scale_color_manual(values = c('Female' = 'Deeppink', 'Male' = 'Deepskyblue')) + geom_abline(intercept = 1, slope = 0) +
  stat_compare_means(method = "wilcox" )
dev.off()

#replot by sex

pdf(paste(Sys.Date(), "ncfem1_to_ppib_nuclear_signal_by_sex.pdf", sep = "_"), height = 5, width = 5)
aggregate %>% ggplot(aes(x=sex, y=value, fill = variable))+ 
  geom_boxplot()+
  scale_fill_manual(values = c("dark red",  "yellow")) +
  geom_point(aes(color = sex),
             position=position_jitterdodge(jitter.width=0.2),
             show.legend = T) +
  scale_color_manual(values = c('Female' = 'Deeppink', 'Male' = 'Deepskyblue')) + geom_abline(intercept = 1, slope = 0) +
  stat_compare_means(method = "wilcox" ) + theme_bw()
dev.off()

##########################3
sink(file = paste(Sys.Date(),"_RT_qPCR_session_Info.txt", sep =""))
sessionInfo()
sink()
