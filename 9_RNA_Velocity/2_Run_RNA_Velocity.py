import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

if __name__ == '__main__':
	
	
	#import sample information
	
	sample_obs = pd.read_csv("cellID_obs.csv")
	umap = pd.read_csv("cell_embeddings.csv")
	
	
	#import loom files
	
	FB1 = anndata.read_loom("KFemale_Blood_FishTEDB_NR.loom")
	FK1 = anndata.read_loom("KFemale_Kidney_FishTEDB_NR.loom")
	FL1 = anndata.read_loom("KFemale_Liver_FishTEDB_NR.loom")
	FS1 = anndata.read_loom("KFemale_Spleen_FishTEDB_NR.loom")
	MB1 = anndata.read_loom("K_Male_Blood_FishTEDB_NR.loom")
	MK1 = anndata.read_loom("K_Male_Kidney_FishTEDB_NR.loom")
	ML1 = anndata.read_loom("K_Male_Liver_FishTEDB_NR.loom")
	MS1 = anndata.read_loom("K_Male_Spleen_FishTEDB_NR.loom")
	
	FB2 = anndata.read_loom("Female_Blood_Cohort2_FishTEDB_NR.loom")
	FK2 = anndata.read_loom("Female_Kidney_Cohort2_FishTEDB_NR.loom")
	FL2 = anndata.read_loom("Female_Liver_Cohort2_FishTEDB_NR.loom")
	MB2 = anndata.read_loom("Male_Blood_Cohort2_FishTEDB_NR.loom")
	ML2 = anndata.read_loom("Male_Liver_Cohort2_FishTEDB_NR.loom")
	MS2 = anndata.read_loom("Male_Spleen_Cohort2_FishTEDB_NR.loom")
	
	FK3 = anndata.read_loom("Female_Kidney_3_FishTEDB_NR.loom")
	FL3 = anndata.read_loom("Female_Liver_3_FishTEDB_NR.loom")
	FS3 = anndata.read_loom("Female_Spleen_3_FishTEDB_NR.loom")
	MK3 = anndata.read_loom("Male_Kidney_3_FishTEDB_NR.loom")
	ML3 = anndata.read_loom("Male_Liver_3_FishTEDB_NR.loom")
	MS3 = anndata.read_loom("Male_Spleen_3_FishTEDB_NR.loom")
	
	
	#rename batches such that the naming scheme is congruent
	
	FB1.obs.index = FB1.obs.index.str.replace("_Fish", "_1_Fish")
	FK1.obs.index = FK1.obs.index.str.replace("_Fish", "_1_Fish")
	FL1.obs.index = FL1.obs.index.str.replace("_Fish", "_1_Fish")
	FS1.obs.index = FS1.obs.index.str.replace("_Fish", "_1_Fish")
	MB1.obs.index = MB1.obs.index.str.replace("_Fish", "_1_Fish")
	MK1.obs.index = MK1.obs.index.str.replace("_Fish", "_1_Fish")
	ML1.obs.index = ML1.obs.index.str.replace("_Fish", "_1_Fish")
	MS1.obs.index = MS1.obs.index.str.replace("_Fish", "_1_Fish")

	FB1.obs.index = FB1.obs.index.str.replace("KFemale", "Female")
	FK1.obs.index = FK1.obs.index.str.replace("KFemale", "Female")
	FL1.obs.index = FL1.obs.index.str.replace("KFemale", "Female")
	FS1.obs.index = FS1.obs.index.str.replace("KFemale", "Female")
	MB1.obs.index = MB1.obs.index.str.replace("K_Male", "Male")
	MK1.obs.index = MK1.obs.index.str.replace("K_Male", "Male")
	ML1.obs.index = ML1.obs.index.str.replace("K_Male", "Male")
	MS1.obs.index = MS1.obs.index.str.replace("K_Male", "Male")
	
	FB2.obs.index = FB2.obs.index.str.replace("Cohort", "")
	FK2.obs.index = FK2.obs.index.str.replace("Cohort", "")
	FL2.obs.index = FL2.obs.index.str.replace("Cohort", "")
	MB2.obs.index = MB2.obs.index.str.replace("Cohort", "")
	ML2.obs.index = ML2.obs.index.str.replace("Cohort", "")
	MS2.obs.index = MS2.obs.index.str.replace("Cohort", "")
	
	#make sure names are contained in objects
	
	sample_obs_FB1 = sample_obs[sample_obs["x"].str.contains("Female_Blood_1_FishTEDB_NR:")]
	sample_obs_FK1 = sample_obs[sample_obs["x"].str.contains("Female_Kidney_1_FishTEDB_NR:")]
	sample_obs_FL1 = sample_obs[sample_obs["x"].str.contains("Female_Liver_1_FishTEDB_NR:")]
	sample_obs_FS1 = sample_obs[sample_obs["x"].str.contains("Female_Spleen_1_FishTEDB_NR:")]
	sample_obs_MB1 = sample_obs[sample_obs["x"].str.contains("Male_Blood_1_FishTEDB_NR:")]
	sample_obs_MK1 = sample_obs[sample_obs["x"].str.contains("Male_Kidney_1_FishTEDB_NR:")]
	sample_obs_ML1 = sample_obs[sample_obs["x"].str.contains("Male_Liver_1_FishTEDB_NR:")]
	sample_obs_MS1 = sample_obs[sample_obs["x"].str.contains("Male_Spleen_1_FishTEDB_NR:")]
	
	sample_obs_FB2 = sample_obs[sample_obs["x"].str.contains("Female_Blood_2_FishTEDB_NR:")]
	sample_obs_FK2 = sample_obs[sample_obs["x"].str.contains("Female_Kidney_2_FishTEDB_NR:")]
	sample_obs_FL2 = sample_obs[sample_obs["x"].str.contains("Female_Liver_2_FishTEDB_NR:")]
	sample_obs_MB2 = sample_obs[sample_obs["x"].str.contains("Male_Blood_2_FishTEDB_NR:")]
	sample_obs_ML2 = sample_obs[sample_obs["x"].str.contains("Male_Liver_2_FishTEDB_NR:")]
	sample_obs_MS2 = sample_obs[sample_obs["x"].str.contains("Male_Spleen_2_FishTEDB_NR:")]
	
	sample_obs_FK3 = sample_obs[sample_obs["x"].str.contains("Female_Kidney_3_FishTEDB_NR:")]
	sample_obs_FL3 = sample_obs[sample_obs["x"].str.contains("Female_Liver_3_FishTEDB_NR:")]
	sample_obs_FS3 = sample_obs[sample_obs["x"].str.contains("Female_Spleen_3_FishTEDB_NR:")]
	sample_obs_MK3 = sample_obs[sample_obs["x"].str.contains("Male_Kidney_3_FishTEDB_NR:")]
	sample_obs_ML3 = sample_obs[sample_obs["x"].str.contains("Male_Liver_3_FishTEDB_NR:")]
	sample_obs_MS3 = sample_obs[sample_obs["x"].str.contains("Male_Spleen_3_FishTEDB_NR:")]
	
	#put into array
	samp_obs = [sample_obs_FB1,sample_obs_FK1,sample_obs_FL1,sample_obs_FS1,sample_obs_MB1,sample_obs_MK1,sample_obs_ML1,sample_obs_MS1,sample_obs_FB2,sample_obs_FK2,sample_obs_FL2,sample_obs_MB2,sample_obs_ML2,sample_obs_MS2,sample_obs_FK3,sample_obs_FL3,sample_obs_FS3,sample_obs_MK3,sample_obs_ML3,sample_obs_MS3]
	
	#put objects into array
	loom_files = [FB1,FK1,FL1,FS1,MB1,MK1,ML1,MS1,FB2,FK2,FL2,MB2,ML2,MS2,FK3,FL3,FS3,MK3,ML3,MS3]
	
	#run function to name things individually
	
	for i in range(len(loom_files)):
		loom_files[i] = loom_files[i][np.isin(loom_files[i].obs.index,sample_obs["x"])]
		loom_files[i] = loom_files[i][np.isin(loom_files[i].obs.index, samp_obs[i])]
		loom_files[i].var_names_make_unique()
	
	all = loom_files[0].concatenate(list(loom_files[1:len(loom_files)]))
	all.obs.index = all.obs.index.str.replace("x.+", "x", regex = True)
	
	all_index = pd.DataFrame(all.obs.index)
	all_index = all_index.rename(columns = {'CellID':'Cell ID'})
	
	#add umap and cluster info
	
	umap = umap.rename(columns = {'Unnamed: 0':'Cell ID'})
	umap_ordered = all_index.merge(umap, on = "Cell ID")
	
	celltype_cols={'B_Cell_Progenitors': '#FD0000',
	"B_cells": "#16FC00",
	"Cholangiocytes": "#1C22F5",
	"Endothelial": "#FDCCCB",
	"Erythrocyte_Progenitors": "#FC16DE",
	"Erythrocytes": "#0DD4FB",
	"Fibroblasts": "#E8E300",
	"Hepatocytes": "#0D600D",
	"Hepatocytes_Efferocytosing": "#FC8B1C",
	"HSPCs": "#623B82",
	"Kidney_distal_tubule": "#AD0056",
	"Kidney_prox_tubule": "#0DFECA",
	"Lymphoid_progenitors": "#D273FF",
	"Macrophages": "#733516",
	"Mast_cells": "#F387D2",
	"Multipotent_progenitors": "#CAEDAD",
	"Myeloid_progenitors": "#73949B",
	"Neutrophil_Progenitors": "#7F9FFD",
	"Neutrophils": "#CC9938",
	"NK_T_cells": "#FE007D",
	"NK_T_progenitor_cells": "#DBC3FF",
	"Thrombocytes": "#7AC91C"}
	
	umap_ordered = umap_ordered.iloc[:,1:]
	all.obsm['X_umap'] = umap_ordered.values
	
	clusters = pd.read_csv("cell_type_obs.csv")
	clusters_ordered = all_index.merge(clusters, on = "Cell ID")
	clusters_ordered = clusters_ordered.rename(columns = {'Cluster':'Cell_Type'})
	clusters_ordered = clusters_ordered.rename(columns = {'Cell ID':'CellID'})
	
	scv.pp.filter_and_normalize(all)
	scv.pp.moments(all)
	scv.tl.velocity(all, mode = 'stochastic')
	scv.tl.velocity_graph(all, n_jobs = 16)
	
	all.obs = all.obs.merge(clusters_ordered, on='CellID')
	
	scv.pl.velocity_embedding(all, basis = 'umap', color = 'Cell_Type', palette=celltype_cols, save = 'all_velo_basic2.pdf')
	scv.pl.velocity_embedding_grid(all, basis = 'umap', color = 'Cell_Type', palette=celltype_cols, save = 'all_velo_grid2.pdf')
	scv.pl.velocity_embedding_stream(all, basis='umap', color = 'Cell_Type', palette=celltype_cols, save = 'all_velo_stream2.pdf')
	scv.pl.velocity_embedding_stream(all, basis='umap', color = 'Cell_Type', palette=celltype_cols,legend_loc='none', save = 'all_velo_stream3.pdf')
	
