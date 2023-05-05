########################################################
	   README - 3_Integration_and_cell_annotation
########################################################

# Software versions are indicated in each cognate script as comments

Run the following scripts in this order:

	1. Import preprocessed data from 3 cohorts and annotate
		- 1_analyze_killifish_Tissues_scRNAseq_3_cohorts_v12.R   : master script to import data from the 3 cohorts, normalize, correct batch effects, annotate cell types
		- 2022-10-10_Zebrafish_cell_markers_forAnnot.txt         : markers used for annotation
		
	2. Plot top markers on UMAP
		- 2_Celltype_marker_and_vitellogenin_plots.R             : script to plot top cell type markers