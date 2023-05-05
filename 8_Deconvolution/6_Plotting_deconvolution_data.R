#Run analysis and make plots from the liver deconvolution data

##########################
#set working directory
setwd("/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/CODE/8_Deconvolution")

###################

library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(scales)

###################
#1. Process deconvolution data for plotting

#read in deconvolution cell type proportion results
decon <- read.table("/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/CODE/8_Deconvolution/output/liver_deconv.txt")

#rename
colnames(decon) <- c("AL1","DR1","DR2","AL2","DR3","DR4","AL3","DR5","AL4","DR6","AL5","AL6","DR7","AL7","DR8","AL8")

#transform dataset
decon <- t(decon)
decon <- decon[rowSums(decon) > 0,]
decon <- as.data.frame(decon)
decon$group <- gsub("Aligned.sortedByCoord.out.bam","",rownames(decon))
decon$sex <- c("M","F","M","F","F","M","M","F","F","M","M","M","F","F","M","F")
decon_m <- melt(decon)
decon_m$group <- gsub("\\d+", "", decon_m$group)
decon_m$group <- factor(decon_m$group, levels = c("AL", "DR"))
decon_m$both <- paste0(decon_m$sex, "_",decon_m$group)

#2. Make pie charts of cell proportions
#set colors
celltype.colors <- c("B_Cell_Progenitors"="#FD0000",
                     "B_cells"="#16FC00",
                     "Cholangiocytes"="#1C22F5",
                     "Endothelial"="#FDCCCB",
                     "Erythrocyte_Progenitors"="#FC16DE",
                     "Erythrocytes"="#0DD4FB",
                     "Fibroblasts"="#E8E300",
                     "Hepatocytes"="#0D600D",
                     "Hepatocytes_Efferocytosing"="#FC8B1C",
                     "HSPCs"="#623B82",
                     "Kidney_distal_tubule"="#AD0056",
                     "Kidney_prox_tubule"="#0DFECA",
                     "Lymphoid_progenitors"="#D273FF",
                     "Macrophages"="#733516",
                     "Mast_cells"="#F387D2",
                     "Multipotent_progenitors"="#CAEDAD",
                     "Myeloid_progenitors"="#73949B",
                     "Neutrophil_Progenitors"="#7F9FFD",
                     "Neutrophils"="#CC9938",
                     "NK_T_cells"="#FE007D",
                     "NK_T_progenitor_cells"="#DBC3FF",
                     "Thrombocytes"="#7AC91C")

#function to get mean values of cell type proportion per group
run_get_mean <- function(x){
  empty_df <- data.frame(matrix(NA, nrow = length(names(table(x$variable))), ncol =2))
  for(i in 1:length(names(table(x$variable)))) {
    cell.type <- subset(x, x$variable == names(table(x$variable))[i])
    empty_df[i,1] <- mean(cell.type$value)
    empty_df[i,2] <- names(table(x$variable))[i]
  }
  colnames(empty_df) <- c("proportion", "cell_type")
  return(empty_df)
}

#get names of groups
group_names <- names(table(decon_m$both))

#initialize vector for plotting
list_for_pie <- list()

#plot pie charts of the average proportion of each cell type per group
for(i in 1:length(group_names)) {
  list_for_pie[[i]] <- subset(decon_m, decon_m$both == group_names[i])
  list_for_pie[[i]] <- run_get_mean(list_for_pie[[i]])
  bp<- ggplot(list_for_pie[[i]], aes(x="", y=proportion, fill=cell_type))+
    geom_bar(width = 1, stat = "identity") + theme_bw() +
    scale_fill_manual(values=celltype.colors)
  pie <- bp + coord_polar("y", start=0)
  pdf(paste(Sys.Date(),group_names[i], "cell_type_pie_deconvolution.pdf", sep = "_"))
  print(pie)
  dev.off()
}

#3. make boxplots

#A. Generate sex differences box plot
#run wilcox.test for significance 
stat.test <- compare_means(
  value ~ sex, data = decon_m, group.by = "variable",
  method = "wilcox.test", ref.group = "F", p.adjust.method = "BH", paired = F,
)
stat.test

# A tibble: 18 × 9
#variable                   .y.   group1 group2       p p.adj p.format p.signif method  
#<fct>                      <chr> <chr>  <chr>    <dbl> <dbl> <chr>    <chr>    <chr>   
#  1 B_Cell_Progenitors         value F      M      0.165   0.37  0.1652   ns       Wilcoxon
#2 B_cells                    value F      M      0.00147 0.013 0.0015   **       Wilcoxon
#3 Cholangiocytes             value F      M      0.161   0.37  0.1605   ns       Wilcoxon
#4 Endothelial                value F      M      0.105   0.31  0.1049   ns       Wilcoxon
#5 Erythrocyte_Progenitors    value F      M      0.00457 0.027 0.0046   **       Wilcoxon
#6 Erythrocytes               value F      M      0.442   0.72  0.4418   ns       Wilcoxon
#7 Fibroblasts                value F      M      0.798   0.96  0.7984   ns       Wilcoxon
#8 Hepatocytes                value F      M      0.645   0.83  0.6454   ns       Wilcoxon
#9 Hepatocytes_Efferocytosing value F      M      0.645   0.83  0.6454   ns       Wilcoxon
#10 Macrophages                value F      M      0.00109 0.013 0.0011   **       Wilcoxon
#11 Mast_cells                 value F      M      0.574   0.83  0.5737   ns       Wilcoxon
#12 Multipotent_progenitors    value F      M      0.0562  0.22  0.0562   ns       Wilcoxon
#13 Myeloid_progenitors        value F      M      0.0621  0.22  0.0621   ns       Wilcoxon
#14 Neutrophil_Progenitors     value F      M      1       1     1.0000   ns       Wilcoxon
#15 Neutrophils                value F      M      0.959   1     0.9591   ns       Wilcoxon
#16 NK_T_cells                 value F      M      0.195   0.39  0.1949   ns       Wilcoxon
#17 NK_T_progenitor_cells      value F      M      0.250   0.45  0.2497   ns       Wilcoxon
#18 Thrombocytes               value F      M      0.878   0.99  0.8785   ns       Wilcoxon

#export deconvolution by sex stats
write.csv(stat.test, file = paste(Sys.Date(), "deconvolution_sex_stats.csv", sep = "_"))

#boxplot with different shapes as different dietary groups
bp <- ggplot(decon_m, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=sex)) + theme_bw()
bp <- bp + stat_pvalue_manual(
  stat.test, x = "variable", y.position = .46,
  label = "p.adj",
  position = position_dodge(0.8)
)
bp <- bp + scale_fill_manual(values = c('F' = 'Deeppink', 'M' = 'Deepskyblue'))
bp.full.sex <- bp + geom_point(aes(group = sex, shape = factor(group)), position=position_jitterdodge(), size = 1)

pdf(paste(Sys.Date(), "sex_diff_decon_full.pdf", sep = "_"), height = 4, width= 8 )
print(bp.full.sex)
dev.off()

#B. Generate sex differences box plot only on significantly differential cell types

sig_cell_types <- subset(decon_m, variable %in% c("B_cells","Erythrocyte_Progenitors","Macrophages"))
stat.test.sig <- subset(stat.test, variable %in% c("B_cells","Erythrocyte_Progenitors","Macrophages"))

#boxplot with different shapes as different dietary groups
bp <- ggplot(sig_cell_types, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=sex)) + theme_bw()
bp <- bp + stat_pvalue_manual(
  stat.test.sig, x = "variable", y.position = .04,
  label = "p.adj",
  position = position_dodge(0.8)
)
bp <- bp + scale_fill_manual(values = c('F' = 'Deeppink', 'M' = 'Deepskyblue'))
bp.sig.sex <- bp + geom_point(aes(group = sex, shape = factor(group)), position=position_jitterdodge(), size = 1)

pdf(paste(Sys.Date(), "sex_diff_decon_full.pdf", sep = "_"), height = 4, width= 5)
print(bp.sig.sex)
dev.off()

###
#C. Assess the effects of dietary restiction

#create a stats test to compare groups split by sex using the t.test and bonferroni test
run_stat_test_dr <- function(x) {compare_means(
  value ~ both, data = x, group.by = "variable",
  method = "wilcox.test", p.adjust.method = "BH", paired = F,
)}

#split by sex
female_only.decon <- subset(decon_m, decon_m$sex == "F")
male_only.decon <- subset(decon_m, decon_m$sex == "M")

#run stats tests and combine
fem_stats_dr <- run_stat_test_dr(female_only.decon)
male_stats_dr <- run_stat_test_dr(male_only.decon)
dr.stats <- rbind(fem_stats_dr, male_stats_dr)
dr.stats

#export dietary restriction stats
write.csv(dr.stats, file = paste(Sys.Date(), "dietary_restriction_stats.csv", sep = "_"))

#plot the data 
bp <- ggplot(decon_m, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=both)) + theme_bw() 
bp <- bp + stat_pvalue_manual(
  dr.stats, x = "variable", y.position = .46,
  label = "p.adj",
  position = position_dodge(c(0.8)))
bp <- bp + scale_fill_manual(values = c('F_AL' = 'Deeppink4', 'M_AL' = 'Deepskyblue4', 'F_DR' = 'Deeppink1', 'M_DR' = 'Deepskyblue1'))
bp.dr <- bp + geom_point(aes(group = both, shape = factor(both)), position=position_jitterdodge(jitter.width = 0), size = 1)

pdf(paste(Sys.Date(), "dietary_restriction_decon.pdf", sep = "_"), height = 4, width= 8 )
print(bp.dr)
dev.off()

#######################
sink(file = paste(Sys.Date(),"_Plotting_deconvolution_session_Info.txt", sep =""))
sessionInfo()
sink()
