########################################################
		        README -- 1_Preprocessing
########################################################

# Software versions are indicated in each cognate script as comments

Run the following scripts in this order:

	0. Prepare the CellRanger reference
		- 0_CellRanger_reference_preparation.txt    : master script describing reference preparation
		- 0b_process_repeat_library_for_10x_v3.pl   : accessory script for reference preparation
		- 0a_gtf_scTE_prep_v2.pl                    : accessory script for reference preparation

	1. run CellRanger for each cohort
		- 1a_cell_ranger_count_Fish_Tissues_withTEs_cohort1.sh  : cohort 1
		- 1b_cell_ranger_count_Fish_tissue_withTEs_cohort2.sh   : cohort 2
		- 1c_cell_ranger_count_Fish_tissue_withTEs_cohort3.sh   : cohort 3
	
	2. Process each cohort with Seurat
	   Removes background RNA with SoupX, performs mitochondrial quality control, stores sample origin information, normalizes gene expression data, performs cell cycle regression, finds and removes doublets with DoubletFinder and scds.
	   - 2a_CLEAN_killifish_Tissues_scRNAseq_cohort1_v4.R
	   - 2b_CLEAN_killifish_Tissues_scRNAseq_cohort2_v4.R
	   - 2c_CLEAN_killifish_Tissues_scRNAseq_cohort3_v8.R
	
