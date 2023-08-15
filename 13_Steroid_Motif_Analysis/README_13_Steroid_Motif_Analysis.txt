########################################################
          README - 13_Steroid_Motif_Analysis
########################################################

# Software versions are indicated in each cognate script as comments

Run the following scripts in this order:

	1. Uses HOMER to generate BED files containing the genomic coordinates of ARE and ERE sites in the turquoise killifish genome.
		- 1_call_ARE_ERE_sites.sh
	
	2. Calls ARE and ERE responsive genes and bins them into groups with a TSS 500 bo, 1 kb, 2 kb downstream of a steroid response element.
		- 2_Parse_Killi_Steroid_Motifs.R
		
	3. Runs GSEA to find sex enrichment for steroid responsive genes within the DESeq2 differentail expression dataset.
		- 3_enrich_ERE_ARE_Killi.R
