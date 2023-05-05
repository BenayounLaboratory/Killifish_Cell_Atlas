#Compare Cell Type Proportions in the Killifish Tissue Atlas

##########################
#set working directory
setwd("/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/CODE/4_Cell_proportions_and_basic_cell_type_stats")

#load libraries
library("scProportionTest")
library(ggplot2)

##########################
#load Seurat object
load("/Volumes/Collaboration_Drive/Killi_Atlas_For_Bryan/Cell_Annotation_Analysis/2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData")

#initialize object for scProportionTest
sc_tissue <- sc_utils(killi.combined.filt)

###############
#function to make comparisons between groups using the scProportionTest

plot_comparisons <- function (groups_of_interest, number_of_bootstraps, title, FDR, FC) {
  comparisons <- data.frame()
  for(i in 2:length(groups_of_interest)) {
    perm_test <- permutation_test(
      sc_tissue, cluster_identity = "Annotation_v1",
      sample_1 = groups_of_interest[[1]], sample_2 = groups_of_interest[[i]],
      sample_identity = "Group", n_permutations = number_of_bootstraps
    )
    perm_df <- perm_test@results$permutation
    perm_df$comp <- paste0(groups_of_interest[[1]],"_vs_",groups_of_interest[[i]])
    colnames(perm_df)[c(2,3)] <- c("Control", "Treatment")
    comparisons <- rbind(comparisons, perm_df)
    comparisons <- comparisons[comparisons$Control > 0 & comparisons$Treatment > 0, ]
  }
  
  #set significance thresholds
  FDR_threshold = FDR
  log2FD_threshold = log2(FC)
  
  ## Mark the significant results.
  suppressWarnings(comparisons[, significance := ifelse(
    FDR < FDR_threshold & obs_log2FD > log2FD_threshold,
    paste("Higher in Males"),
    ifelse(FDR < FDR_threshold & obs_log2FD < -log2FD_threshold,
           paste("Higher in Females"),
           "n.s.")
  )
  ]
  )
  
  #plot
  plot <- ggplot(comparisons, aes(x = clusters, y = obs_log2FD)) +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = factor(significance) , shape = comp)) +
    theme_bw() +
    geom_hline(yintercept = log2FD_threshold, lty = 2) +
    geom_hline(yintercept = -log2FD_threshold, lty = 2) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = c("Higher in Males" = "Deepskyblue",
                                  "Higher in Females" = "Deeppink",
                                  "n.s." = "gray")) +
    coord_flip() + ggtitle(title) +
    ylim(-4,4)
  
  return(plot)
}

plot_comparisons_sex <- function (groups_of_interest, number_of_bootstraps, title, FDR, FC) {
  comparisons <- data.frame()
  for(i in 2:length(groups_of_interest)) {
    perm_test <- permutation_test(
      sc_tissue, cluster_identity = "Annotation_v1",
      sample_1 = groups_of_interest[[1]], sample_2 = groups_of_interest[[i]],
      sample_identity = "Sex", n_permutations = number_of_bootstraps
    )
    perm_df <- perm_test@results$permutation
    perm_df$comp <- paste0(groups_of_interest[[1]],"_vs_",groups_of_interest[[i]])
    colnames(perm_df)[c(2,3)] <- c("Control", "Treatment")
    comparisons <- rbind(comparisons, perm_df)
    comparisons <- comparisons[comparisons$Control > 0 & comparisons$Treatment > 0, ]
  }
  
  #set significance thresholds
  FDR_threshold = FDR
  log2FD_threshold = log2(FC)
  
  ## Mark the significant results.
  suppressWarnings(comparisons[, significance := ifelse(
    FDR < FDR_threshold & obs_log2FD > log2FD_threshold,
      paste("Higher in Males"),
        ifelse(FDR < FDR_threshold & obs_log2FD < -log2FD_threshold,
          paste("Higher in Females"),
          "n.s.")
  )
  ]
  )
  
  #plot
  plot <- ggplot(comparisons, aes(x = clusters, y = obs_log2FD)) +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = factor(significance) , shape = comp)) +
    theme_bw() +
    geom_hline(yintercept = log2FD_threshold, lty = 2) +
    geom_hline(yintercept = -log2FD_threshold, lty = 2) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = c("Higher in Males" = "Deepskyblue",
                                  "Higher in Females" = "Deeppink",
                                  "n.s." = "gray")) +
    coord_flip() + ggtitle(title) +
    ylim(-4,4)
  
  return(plot)
}
###########################
#generate a list of comparisons. These should be taken from the metadata variable "
#crucially- keep the first item as the baseline comparison i.e. if comparing to a young baseline, add the young variable first

#by organ
blood <-list("Blood_F", "Blood_M")
liver <-list("Liver_F", "Liver_M")
kidney <-list("Kidney_F", "Kidney_M")
spleen <-list("Spleen_F", "Spleen_M")

#by sex
sex <- list("F", "M")

#################
#run the plots

#input the list, bootstraps, plot title, significance threshold value (p val), fold change threshold value

pdf(paste0(Sys.Date(),"_GRZ_Blood_Cell_Proportions.pdf"))
plot_comparisons(blood,  100, "Female Baseline Blood Cell Proportions", 0.05, 2)
dev.off()

pdf(paste0(Sys.Date(),"_GRZ_Liver_Cell_Proportions.pdf"))
plot_comparisons(liver,  100, "Female Baseline Liver Cell Proportions", 0.05, 2)
dev.off()

pdf(paste0(Sys.Date(),"_GRZ_Kidney_Cell_Proportions.pdf"))
plot_comparisons(kidney, 100, "Female Baseline Kidney Cell Proportions", 0.05, 2)
dev.off()

pdf(paste0(Sys.Date(),"_GRZ_Spleen_Cell_Proportions.pdf"))
plot_comparisons(spleen, 100, "Female Baseline Spleen Cell Proportions", 0.05, 2)
dev.off()

#input the list, bootstraps, plot title, significance threshold value (p val), fold change threshold value for all tissues
pdf(paste0(Sys.Date(),"_GRZ_All_Tissues_Cell_Proportions.pdf"))
plot_comparisons_sex(sex, 100, "Female Baseline Blood Atlas Cell Proportions", 0.05, 2)
dev.off()

#######################
sink(file = paste(Sys.Date(),"_Celltype_Proportions_session_Info.txt", sep =""))
sessionInfo()
sink()
