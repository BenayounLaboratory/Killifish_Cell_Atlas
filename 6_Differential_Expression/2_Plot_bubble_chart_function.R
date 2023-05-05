library(ggplot2) 
library(scales) 
theme_set(theme_bw())


# my.gsea.data <- go.results[[8]]
# max.path.plot = 20

get_sex_bubble_plot <- function(my.cell.type, my.gsea.data, max.path.plot = 20){
  
  my.gsea.pos <- my.gsea.data[my.gsea.data$NES > 0,]
  my.gsea.neg <- my.gsea.data[my.gsea.data$NES < 0,]
  
  my.pos.sort <- sort(my.gsea.pos$NES, index.return = T, decreasing = T) # largest value is top (positive)
  my.neg.sort <- sort(my.gsea.neg$NES, index.return = T, decreasing = F) # largest value is top (negative)
  
  if ( (nrow(my.gsea.pos) > round(max.path.plot/2)) && (nrow(my.gsea.neg) > round(max.path.plot/2)) ) {
    
    # if enough on both sides
    my.gsea.sex.2 <- rbind(my.gsea.pos[my.pos.sort$ix[1:round(max.path.plot/2)],],
                           my.gsea.neg[my.neg.sort$ix[1:round(max.path.plot/2)],])
    
  } else {
    
    # if not enough on both sides
    my.gsea.sex.2 <- rbind(my.gsea.pos[my.pos.sort$ix[1:(min(round(max.path.plot/2),nrow(my.gsea.pos)))],],
                           my.gsea.neg[my.neg.sort$ix[1:(min(round(max.path.plot/2),nrow(my.gsea.neg)))],])
    
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
  
  # Female/Male color scale
  my.max <-  3 # add a value in case there are not terms biased in one sex
  my.min <- -3 # add a value in case there are not terms biased in one sex
  
  my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
  my.scaled <- rescale(my.values, to = c(0, 1))
  my.color.vector.sex <- c("deepskyblue","lightskyblue","lightskyblue1","lightcyan","white","lavenderblush","plum1","orchid1","deeppink")
  
  my.plot <- ggplot(my.gsea.sex.sorted,aes(x=CellType,y=PathName,colour=NES,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
  my.plot <- my.plot + ggtitle("GOBP gsea Analysis") + labs(x = "-log10(pvalue)", y = "")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.sex, 
                                              na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(-3, 3) )
  my.plot
  
  my.pdfname <- paste(Sys.Date(),"GSEA_SEX_BALLOON_plot",my.cell.type,"top", nrow(my.gsea.sex.sorted),"significant_GOBP_terms.pdf", sep="_")
  
  pdf(my.pdfname, onefile=T, height = 6, width= 10 )
  print(my.plot)
  dev.off()  
  
  
}

