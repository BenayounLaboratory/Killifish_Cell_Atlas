########################################################
		        README -- 8_Deconvolution
########################################################

# Software versions are indicated in each cognate script as comments

Run the following scripts in this order:

	1. Trims adapters from reads downloaded from SRA.
		- 1_Trim_adapters.sh
	
	2. Aligns reads to the turquoise killifish genome using STAR.
		- 2_STAR_align.sh
		
	3. Counts aligned reads using featureCounts.
		- 3_Count_reads.sh

	4. Retains only cells from the liver from the tissue atlas Seurat object.
		- 4_Liver_isolation_for_deconvolution.R

	5. Deconvolutes liver bulk RNA-seq data using the liver cells from the tissue atlas using CSCDRNA. Outputs a table of cell proportions for each bulk RNA-seq replciate.
		- 5_Liver_deconvolution.R
		
	6. Generates pie charts of cell proportions for each group from McKay et al., boxplots of cell proportions by sex and dietary group.
		- 6_Plotting_deconvolution_data.R