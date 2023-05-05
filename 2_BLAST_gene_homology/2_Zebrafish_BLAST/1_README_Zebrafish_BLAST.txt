# 2022-02-11

# get protein sequences from NCBI for killifish 2015 assembly
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/465/895/GCF_001465895.1_Nfu_20140520/

# get Zebrafish proteins from Ensembl Biomart 105

		Dataset 30368 / 37241 Genes
		Zebrafish genes (GRCz11)

		Filters
			Gene type: IG_C_gene , IG_C_pseudogene , IG_J_pseudogene , IG_pseudogene , IG_V_pseudogene , processed_pseudogene , protein_coding , pseudogene , TR_D_gene , TR_J_gene , TR_V_gene
		
		Attributes
			Peptide
			Gene name
			Gene stable ID
			Protein stable ID
			Transcript stable ID
		
				2022-02-11_Zebrafish_Proteome_Ens105.fa


BLASTn tabular output format 6
Column headers: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore



# reverse blast
makeblastdb -in GCF_001465895.1_Nfu_20140520_protein.faa -out Nfu_20140520_protein -dbtype prot -title "Killi 2015 Database"
blastp -query 2022-02-11_Zebrafish_Proteome_Ens105.fa  -db Nfu_20140520_protein -evalue 1e-5 -out 2022-02-22_Zebrafish_Nfur2015_BestHits_1e-5_REVERSE.txt -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -outfmt 6 -max_target_seqs 1
