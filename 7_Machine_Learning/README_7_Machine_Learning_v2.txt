########################################################
              README - 7_Machine_Learning
########################################################

# Software versions are indicated in each cognate script as comments

Run the following scripts in this order:

	1. Contains functions to run balanced accuracy tests on machine-learning models.
		- 1_ML_helper_functions.R

	2. Trains, tests, and validates RF and GBM machine-learning models on tissue atlas to predict sex from cell-level gene expression data. 
	   Uses functions from 1_ML_helper_functions.R
	   Trains models for each cell type with enough cells in Cohort 1 (> 100 F + 100 M)
		- 2_run_sex_ML_killi_cell_atlas_PER_CELL_TYPE.R
