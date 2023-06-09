########################################################
          README - 6_Differential_Expression
########################################################

# Software versions are indicated in each cognate script as comments

Run the following scripts in this order:

	1. Function to generate Jaccard plots used in 3_analyze_DE_genes_MUSCAT_Pseudobulk_NEW_annot.R
		- 1_Plot_Jaccard_function.R
	
	2. Function to generate GO bubble plots used in 3_analyze_DE_genes_MUSCAT_Pseudobulk_NEW_annot.R
		- 2_Plot_bubble_chart_function.R
		
	3. Runs muscat to pseudobulk cell types, RUV to reduce batch/tissue effects on gene expression, runs GO BP analysis with ClusterProfiler, and generates Jaccard index plots.
		- 3_analyze_DE_genes_MUSCAT_Pseudobulk_NEW_annot.R
		
	4. Plots GO results for each cell type with standardized NES and FDR scales.
		- 4_Cell_type_GO_Plotting.R