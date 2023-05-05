#Get the median unique genes per cell type per sex and generate pie charts of cell type proportions by sex and tissue

##########################
#set working directory
setwd("/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/CODE/4_Cell_proportions_and_basic_cell_type_stats")

# Load libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggbeeswarm)

# get UMAP
load("/Volumes/Collaboration_Drive/Killi_Atlas_For_Bryan/Cell_Annotation_Analysis/2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData")

#1. Plot cell type proportions as pie chart
#get sexes and tissues
sexes <- unique(killi.combined.filt@meta.data$Sex)
tissues <- unique(killi.combined.filt@meta.data$Tissue)

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

#generate pie charts of cell type proportions by sex and tissue
plotter <- function(a,b){
  sub <- subset(killi.combined.filt, subset = Sex %in% a & Tissue %in% b)
  cell.comp <- as.data.frame(table(sub@meta.data$Annotation_v1))
  cell.comp <- subset(cell.comp, cell.comp$Freq > 10)
  return(cell.comp)
}
empty <- list()
i <- 0
while(i <= (length(sexes) * length(tissues))){
  for(x in 1:length(sexes)) {
    for(y in 1:length(tissues)){
      output <- plotter(sexes[[x]], tissues[[y]])
      colnames(output)[1] <- c(paste0(sexes[[x]], "_", tissues[[y]]))
      i <- i+1
      empty[[i]] <-output
    }
  }
}
plotlist <- list()
for(j in 1:(length(empty)/2)){
  plotlist[[j]]  <- empty[[j]] %>% 
    arrange(desc(empty[[j]][1])) %>%
    mutate(prop = Freq / sum(empty[[j]]$Freq) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
}
fun_plot <- function(k) {
  ggplot(k, aes(x="", y=prop, fill=k[,1])) +
    geom_bar(stat="identity") +
    coord_polar("y", start=0) +
    theme_void() + ggtitle(colnames(k)[1]) +
    theme(legend.position="none") + scale_fill_manual(values=celltype.colors) +
    geom_text(aes(y = ypos, label = k[,1]), color = "black", size=3)
}

#print pdf of pie chart
pdf(paste(Sys.Date(),"cell_type_pie.pdf", sep = "_"))
ggarrange(plotlist = lapply(plotlist, fun_plot), nrow = 2, ncol = 4)
dev.off()

#2. Plot median terms per cell type per batch per sex per cohort
sexes <- unique(killi.combined.filt@meta.data$Sex)
batches <- unique(killi.combined.filt@meta.data$Batch)

gene.exp <- NULL
for(x in 1:length(sexes)) {
  for(y in 1:length(batches)) {
    sub.sb <- subset(killi.combined.filt, subset = Sex %in% sexes[[x]] & Batch %in% batches[[y]])
    ct.names <- names(table(sub.sb@meta.data$Annotation_v1))
    for(z in 1:length(ct.names)) {
      sub.ct <- subset(sub.sb, subset = Annotation_v1 %in% ct.names[[z]])
      data <- median(sub.ct@meta.data[["nFeature_RNA"]])
      data <- as.data.frame(data)
      data[,2] <- sexes[[x]]
      data[,3] <- batches[[y]]
      data[,4] <- ct.names[[z]]
      colnames(data) <- c("median_nFeatures", "sex", "batch", "celltype")
      gene.exp <- rbind(gene.exp, data)
    }
  }
}

#plot median number of genes expressed per cell type per cohort
pdf(paste(Sys.Date(),"median_genes_per_cell_type.pdf", sep = "_"))
ggplot(gene.exp, aes(x=median_nFeatures, y=celltype, fill=sex)) +
  geom_boxplot() + theme_bw() + geom_beeswarm(dodge.width = 0.75) +
  scale_fill_manual(values=c("DeepPink", "DeepSkyBlue"))
dev.off()

#######################
sink(file = paste(Sys.Date(),"_Plotting_cell_types_session_Info.txt", sep =""))
sessionInfo()
sink()

