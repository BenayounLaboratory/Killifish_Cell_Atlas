# logFC.data = my.logFC.data
# FDR.data <- my.FDR.data
# fdr.thres = 0.05

calc_jaccard <- function(logFC.data, FDR.data, fdr.thres = 0.05){
  
  my.up  <- logFC.data > 0
  my.sig <- FDR.data   < fdr.thres
  
  up.genes   <- data.frame(matrix(0,length(my.exp.genes),length(deseq.res.list.genes)))
  colnames(up.genes) <- names(deseq.res.list.genes)
  rownames(up.genes) <- my.exp.genes
  
  down.genes <- data.frame(matrix(0,length(my.exp.genes),length(deseq.res.list.genes)))
  colnames(down.genes) <- names(deseq.res.list.genes)
  rownames(down.genes) <- my.exp.genes
  
  for (i in 1:length(deseq.res.list.genes)) {
    for (j in 1:length(my.exp.genes)) {
      
      up.genes[j,i]   <- bitAnd(my.up[j,i], my.sig[j,i])
      down.genes[j,i] <- bitAnd(!my.up[j,i], my.sig[j,i])
      
    }
  }
  up.genes[is.na(up.genes)]     <- 0
  down.genes[is.na(down.genes)] <- 0
  
  my.jaccard.up <- matrix(0,ncol(up.genes),ncol(up.genes))
  rownames(my.jaccard.up) <- colnames(up.genes)
  colnames(my.jaccard.up) <- colnames(up.genes)
  
  for (i in 1:ncol(up.genes) ) {
    for (j in 1:ncol(up.genes) ) {
      
      my.inter <- sum(bitAnd(up.genes[,i],up.genes[,j]))
      my.union <- sum(bitOr(up.genes[,i],up.genes[,j]))
      
      my.jaccard.up[i,j] <- my.inter/my.union
    }
  }
  
  
  my.jaccard.dwn <- matrix(0,ncol(down.genes),ncol(down.genes))
  rownames(my.jaccard.dwn) <- colnames(down.genes)
  colnames(my.jaccard.dwn) <- colnames(down.genes)
  
  for (i in 1:ncol(down.genes) ) {
    for (j in 1:ncol(down.genes) ) {
      
      my.inter <- sum(bitAnd(down.genes[,i],down.genes[,j]))
      my.union <- sum(bitOr(down.genes[,i],down.genes[,j]))
      
      my.jaccard.dwn[i,j] <- my.inter/my.union
    }
  }
  
  return(list(my.jaccard.up, my.jaccard.dwn))
}

