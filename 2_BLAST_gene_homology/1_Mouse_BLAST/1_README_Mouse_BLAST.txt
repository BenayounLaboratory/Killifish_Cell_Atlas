# 2022-10-10

# get protein sequences from NCBI for killifish 2015 assembly
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/465/895/GCF_001465895.1_Nfu_20140520/

# get Mouse proteins from Ensembl Biomart 107

		Dataset 22464 / 56748 Genes
		Mouse genes (GRCm39)

		Filters
			Gene type: Gene type: IG_C_gene , IG_D_gene , IG_J_gene , IG_LV_gene , IG_V_gene , protein_coding , TR_C_gene , TR_D_gene , TR_J_gene , TR_V_gene
		
		Attributes
			Peptide
			Gene name
			Gene stable ID
			Protein stable ID
			Transcript stable ID
		
				2022-10-10_Mouse_Proteome_Ens107.fa


BLASTn tabular output format 6
Column headers: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

# reverse blast to gget mouse homologs
makeblastdb -in GCF_001465895.1_Nfu_20140520_protein.faa -out Nfu_20140520_protein -dbtype prot -title "Killi 2015 Database"
blastp -query 2022-10-10_Mouse_Proteome_Ens107.fa  -db Nfu_20140520_protein -evalue 1e-3 -out 2022-10-10_Mouse_Nfur2015_BestHits_1e-3_REVERSE.txt -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -outfmt 6 -max_target_seqs 1
