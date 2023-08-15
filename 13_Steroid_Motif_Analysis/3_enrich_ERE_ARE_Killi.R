setwd('/Users/berenice/Dropbox/Manuscripts_and_Publications/2023/2023_Killifish_single_cell_Atlas_Bryan/TEST/Steroid_motif')
options(stringsAsFactors = F)

library(clusterProfiler)
library(readxl)
library(ggplot2)
library(scales) 

# 2023-06-08

########################################
# 1. read ARE/ERE killi gene set (HOMER calls)
ere_are.gset <- read.gmt('2023-06-08_HOMER_2015_Killi_Genome_ERE_ARE_sites_Promoters.gmt')

########################################
# 2. Read DEseq2 results from Table S4
names.deseq.res <- data.frame(read_xlsx('TableS4_DESeq2_results_by_cell_type_v1.xlsx', sheet = "Index", skip = 2, col_names = F))

my.cell.names <- unlist(lapply(strsplit(names.deseq.res[,2], " from ", fixed = T),'[[',2))
my.cell.names <- gsub("/","_",my.cell.names)
my.cell.names <- gsub("-",".",my.cell.names)
my.cell.names <- gsub(" ","_",my.cell.names)

# Create list object to read DEseq2 results
deseq2.res <- vector(mode = "list", length = length(my.cell.names))
names(deseq2.res) <- my.cell.names

for (i in 1:length(deseq2.res)) {
  deseq2.res[[i]] <- data.frame(read_xlsx('TableS4_DESeq2_results_by_cell_type_v1.xlsx', sheet = i + 1))
}

########################################
# 3.  Loop over Deseq2 results and run GSEA

# Create list object to receive GSEA results
gsea.results <- vector(mode = "list", length = length(deseq2.res))
names(gsea.results) <- names(deseq2.res)

cprof.results <- vector(mode = "list", length = length(deseq2.res))
names(cprof.results) <- names(deseq2.res)

set.seed(1234567)

for (i in 1:length(deseq2.res)) {
  
  # get outprefix
  my.outprefix <- paste0(Sys.Date(),"_ClusterProfiler_GSEA_ARE_ERE",names(deseq2.res)[i])
  
  # Prepare GeneList using log2FC to rsank genes
  sc.glist        = deseq2.res[[i]]$log2FoldChange
  names(sc.glist) = deseq2.res[[i]]$...1
  sc.glist        = sort(sc.glist, decreasing = TRUE)
  
  # GSEA
  ere_are.gsea <- GSEA(geneList     = sc.glist,
                       TERM2GENE    = ere_are.gset,
                       minGSSize    = 5 ,
                       maxGSSize    = 5000,
                       verbose      = FALSE,
                       pvalueCutoff = 1,
                       seed = T)
  View(ere_are.gsea@result)
  
  # write results ere_are.gsea file
  write.table(ere_are.gsea@result, file = paste0(my.outprefix,"_GSEA_Analysis_ALL.txt"), quote = F, sep = "\t")
  
  # save results into object
  gsea.results[[i]]  <- ere_are.gsea@result
  cprof.results[[i]] <- ere_are.gsea
  
}

########################################
# 4. make summary plot

###################### Summary plot across cell types ######################  

my.pathways <- rownames(gsea.results[[i]])


############
# prepare output data
# p-val matrix
my.matrix <- matrix(0,length(my.pathways),length(gsea.results)) # default: -log10(1) pval == 0 no enrichment

# Enrichment matrix
my.matrix2 <- matrix(0,length(my.pathways),length(gsea.results)) # initialize with Enrichment = 0 if no enrich

# matrix with record of significance
my.matrix3 <- matrix(0,length(my.pathways),length(gsea.results)) # to get sigificant pathways

colnames(my.matrix)  <- names(gsea.results)
colnames(my.matrix2) <- names(gsea.results)
colnames(my.matrix3) <- names(gsea.results)
rownames(my.matrix)  <- my.pathways
rownames(my.matrix2) <- my.pathways
rownames(my.matrix3) <- my.pathways

# collect data from clusterprofiler run
for (i in 1:length(my.pathways)) {
  for (j in 1:length(gsea.results)) { # tissues 
    # determine position
    my.id <- which(gsea.results[[j]]$ID %in% my.pathways[i])
    
    if(length(my.id) == 1) { # if was present in this tissue
      # extract stats: FDR, NES and presence
      my.matrix[i,j] <- -log10(gsea.results[[j]]$p.adjust[my.id]) # log(0) is undefined
      my.matrix2[i,j] <- gsea.results[[j]]$NES[my.id]
      my.matrix3[i,j] <- 1
    }
  }
}


# Plot ALL
my.res.enrich   <- data.frame(my.matrix2)
my.pval.enrich  <- data.frame(my.matrix )
my.res.enrich$Pathnames <- rownames(my.res.enrich)

# format for ggplot
my.res.enrich3 <- cbind(my.res.enrich[,c('Pathnames',names(gsea.results)[1])],rep(names(gsea.results)[1],dim(my.res.enrich)[1]),my.pval.enrich[,names(gsea.results)[1]])
colnames(my.res.enrich3) <- c('Pathnames','NES','condition','minusLog10Pval')
for ( h in 2:length(names(gsea.results))) {
  my.new <- cbind(my.res.enrich[,c('Pathnames',names(gsea.results)[h])],rep(names(gsea.results)[h],dim(my.res.enrich)[1]),my.pval.enrich[,names(gsea.results)[h]])
  colnames(my.new) <- colnames(my.res.enrich3)
  my.res.enrich3 <- rbind(my.res.enrich3, 
                          my.new)
}

# Female/Male color scale
my.max <-  2.5 # add a value in case there are not terms biased in one sex
my.min <- -2.5 # add a value in case there are not terms biased in one sex
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.color.vector.sex <- c("deepskyblue","lightskyblue","lightskyblue1","lightcyan","white","lavenderblush","plum1","orchid1","deeppink")

# to preserve the wanted order
my.res.enrich3$condition <- factor(my.res.enrich3$condition, levels = unique(my.res.enrich3$condition))
my.res.enrich3$Pathnames <- factor(my.res.enrich3$Pathnames, levels = rev(unique(my.res.enrich3$Pathnames)))

my.pdfname <- paste(Sys.Date(),"Enrichment_BALLOON_plot_scRNAseq_killi_SexDim_ARE_ERE_enrichment_analysis.pdf", sep="_")

pdf(my.pdfname, onefile=F, height = 5, width=8)
my.plot <- ggplot(my.res.enrich3,aes(x=condition,y=Pathnames,colour=NES,size=minusLog10Pval))+ theme_bw()+ geom_point(shape = 16) 
my.plot <- my.plot + ggtitle("Estrogen and Androgen response elements") + labs(x = "Tissue/condition", y = "Gene Set")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.sex,na.value = "grey50", guide = "colourbar", values = my.scaled, limits = c(-2.5, 2.5) )
my.plot <- my.plot + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
my.plot <- my.plot + scale_size_area(breaks = c(0,2,4) )
print(my.plot)
dev.off()  

