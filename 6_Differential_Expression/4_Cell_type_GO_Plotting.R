#setwd
setwd("/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/CODE/5_Differential_Expression")

#import libraries
library(ggplot2)
library(ggpubr)
library(scales)

#import GO Results
load("/Volumes/Collaboration_Drive/Killi_Atlas_For_Bryan/PseudoBulk_Differential_Expression/2023-02-10_pseudobulk_killi_cell_types_SEX_DIM_GOBP_GSEA.RData")

#write function in which the top 5 GO terms are taken for each cell type
get_sex_bubble_plot <- function(my.cell.type, my.gsea.data, max.path.plot = 20){
  my.gsea.pos <- my.gsea.data[my.gsea.data$NES > 0,]
  my.gsea.neg <- my.gsea.data[my.gsea.data$NES < 0,]
  my.pos.sort <- sort(my.gsea.pos$NES, index.return = T, decreasing = T) # largest value is top (positive)
  my.neg.sort <- sort(my.gsea.neg$NES, index.return = T, decreasing = F) # largest value is top (negative)
  if ( (nrow(my.gsea.pos) > round(max.path.plot/4)) && (nrow(my.gsea.neg) > round(max.path.plot/2)) ) {
    # if enough on both sides
    my.gsea.sex.2 <- rbind(my.gsea.pos[my.pos.sort$ix[1:round(max.path.plot/4)],],
                           my.gsea.neg[my.neg.sort$ix[1:round(max.path.plot/4)],])
  } else {
    # if not enough on both sides
    my.gsea.sex.2 <- rbind(my.gsea.pos[my.pos.sort$ix[1:(min(round(max.path.plot/4),nrow(my.gsea.pos)))],],
                           my.gsea.neg[my.neg.sort$ix[1:(min(round(max.path.plot/4),nrow(my.gsea.neg)))],])
    # if only pathways on one side, there is an NI to remove
    my.gsea.sex.2 <- my.gsea.sex.2[!is.na(my.gsea.sex.2$ID),]
  }
  # create -log10 FDR for plotting
  my.gsea.sex.2$minlog10fdr  <- -log10(my.gsea.sex.2$qvalues)
  my.gsea.sex.2$PathName  <- paste0(my.gsea.sex.2$ID, " ", my.gsea.sex.2$Description)
  my.sorting <- sort(my.gsea.sex.2$minlog10fdr, index.return = T, decreasing = T)
  my.gsea.sex.sorted <- my.gsea.sex.2[my.sorting$ix,]
  # create and preserve wanted display order
  my.gsea.sex.sorted$sex <- ifelse(my.gsea.sex.sorted$NES < 0, "Male", "Female")  # male/female avg flag
  my.gsea.sex.sorted <- my.gsea.sex.sorted[order(my.gsea.sex.sorted$sex), ]
  my.max.char <- max(nchar(my.gsea.sex.sorted$Description))
  my.gsea.sex.sorted$PathName <- factor(my.gsea.sex.sorted$PathName, levels = rev(unique(my.gsea.sex.sorted$PathName)))
  my.gsea.sex.sorted$CellType <- factor(rep("PB",length(  my.gsea.sex.sorted$PathName)))
  return(my.gsea.sex.sorted)
}

#initialize go_list list
go_list <- list()

#for loop to add GO result if the cell type passed significance
for(i in 1:length(go.results)) {
  if(nrow(go.results[[i]])>0) { # only if something passed significance
    go_list[[i]] <- get_sex_bubble_plot(names(go.results)[i], go.results[[i]])
    go_list[[i]]$CellType <- names(go.results)[[i]] #append new category to add name of cell type
  }
}

#retaining only the significant results
go_list <- go_list[c(1:3,5:7,10,11,13)]

# Set the Female/Male color scale
my.max <-  3 # add a value in case there are not terms biased in one sex
my.min <- -3 # add a value in case there are not terms biased in one sex
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector.sex <- c("deepskyblue","lightskyblue","lightskyblue1","lightcyan","white","lavenderblush","plum1","orchid1","deeppink")

#function to plot using a consistent NES and log10FDR scale from the min to max log10FDR value (1.4 to 1.86)
plot_celltype_GO <- function(x){
  my.plot <- ggplot(x,aes(x=CellType,y=PathName,colour=NES,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)+scale_size_continuous(limits=c(1.4,1.86))
  my.plot <- my.plot + labs(x = "-log10(pvalue)", y = "")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.sex, 
                                              na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(-3, 3) )  + theme(text = element_text(size = 7))
  my.plot <- my.plot + ggtitle(x[1,15]) #add in title of cell type
  my.plot
}

#apply plot function to list of cell type GO terms
pdf(paste(Sys.Date(),"significant_GO_results_by_cell_type.pdf", sep = "_"), onefile=T, height = 8, width= 15 )
ggarrange(plotlist = lapply(go_list, plot_celltype_GO), nrow = 3, ncol = 3)
dev.off()

#######################
sink(file = paste(Sys.Date(),"_Cell_type_GO_plots_Session_info.txt", sep =""))
sessionInfo()
sink()

