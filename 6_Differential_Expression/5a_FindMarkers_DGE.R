#run FindMarkers to see if certain cell types are more differential

setwd("/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/REVISION_ANALYSES/Findmarkers_DGE")

#load libraries
library(Seurat)

#load the tissue atlas
load("/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/REVISION_ANALYSES/edgeR_LRT/2023-02-10_Killifish_Tissues_3_cohorts10XGenomics_Singlets_Seurat_object_With_Annotation_Labels_QC_CellsONLY.RData")

#get the names of the cell types
names_of_cells <- names(table(killi.combined.filt@meta.data[["Annotation_v1"]]))

#run FindMarkers on each cell type and store in a list
DGE_list <- list()
for(i in 1:length(names_of_cells)) {
  diff_list[[i]] <- FindMarkers(subset(killi.combined.filt, subset = Annotation_v1 %in% names_of_cells[i]),
                                ident.1 = "F", ident.2 = "M", group.by = "Sex", min.pct = -Inf,
                                logfc.threshold = -Inf, min.cells.feature = 1, min.cells.group = 1)}

#set names of DGE list items
names(DGE_list) <- names_of_cells

#write results to table
for(i in 1:length(names_of_cells)) {
write.table(DGE_list[[i]], file = paste(Sys.Date(), names(DGE_list)[i], "FindMarkers_DGE.txt", sep = "_"), sep = "\t" , row.names = T, quote = F)
}

#keep only cell types that are differently expressed per DESeq2 for comparison
DGE_list <- DGE_list[c(2,4,5,6,8,9,14,15,16,17,19,20,22)]

# get color palette used for UMAP
colors.cells <- read.table('/Users/bryanteefy/Dropbox/2022_Killifish_cell_Atlas/REVISION_ANALYSES/edgeR_LRT/2023-04-12_color_palette_annotation.txt', header = T)

# get correspondence cell type/color
ct.unique <- data.frame("label" = colors.cells$cell_type, "cols" = colors.cells$col)

# filter to EdgeR tested cell types only
ct.unique           <- ct.unique[ct.unique$label %in% names(DGE_list),]
rownames(ct.unique) <- ct.unique$label

## Order by pvalue:
sex.results <- lapply(DGE_list,function(x) {x[order(x$p_val_adj),]})
n           <- sapply(sex.results, nrow)
names(n)    <- names(sex.results)

cols <- list()
xlab <- character(length = length(sex.results))
for(i in seq(along = sex.results)){
  cols[[i]] <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i <- sex.results[[i]]$p_val_adj < 0.05
  cols[[i]][ind.sig.i] <- ct.unique[names(sex.results)[i], "cols"]
  xlab[i] <- paste(names(sex.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(sex.results)

#plot
pdf(paste0(Sys.Date(),"_stripplot_FindMarkers_with_cell_type_colors.pdf"), width = 6, height = 5)
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 13.5),
     ylim = c(-3, 3),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change (F / M)"
)
abline(h = 0)
abline(h = seq(-5, 5, by = 1)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(sex.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(sex.results[[i]])), amount = 0.2),
         y = rev(sex.results[[i]]$avg_log2FC),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]),
         cex = 0.7)
}
axis(1,
     at = 1:13,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-5, 5, 1))
box()
dev.off()

#######################
sink(file = paste(Sys.Date(),"FindMarkers_Per_Cell_Type.txt", sep =""))
sessionInfo()
sink()


