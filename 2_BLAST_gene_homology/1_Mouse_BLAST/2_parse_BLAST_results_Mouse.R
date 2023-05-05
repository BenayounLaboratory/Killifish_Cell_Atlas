setwd('/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2022-02-01_Killifish_scRNAseq_ATLAS_analyses/Mouse_alignment/')
options(stringsAsFactors = FALSE)

# 2022-10-11
# Parse Killifish/Mouse BLAST alignments for homology mapping

#########################################################################################################################
# read in result table and add in column names
blast.rev.res <- read.table('2022-10-10_Mouse_Nfur2015_BestHits_1e-3_REVERSE.txt', sep = "\t", header = F)
colnames(blast.rev.res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast.rev.res$Mmu_Symbol  <- unlist(lapply(strsplit(blast.rev.res$qseqid,"|", fixed = TRUE),'[',1))
blast.rev.res$Mmu_ENSMUSP <- unlist(lapply(strsplit(blast.rev.res$qseqid,"|", fixed = TRUE),'[',2))
blast.rev.res$Mmu_ENSMUST <- unlist(lapply(strsplit(blast.rev.res$qseqid,"|", fixed = TRUE),'[',3))

# read in killifish gene annotation table
killi.annot    <- read.csv('GCF_001465895.1_Nfu_20140520_feature_table.txt', header = F, sep = "\t")
killi.annot.cl <- killi.annot[,c(1, 11, 14, 15)]
colnames(killi.annot.cl) <- c("Type", "Nfur_Accession", "Nfur_Description", "Nfur_Symbol")

# select on protein only
killi.annot.cl <- killi.annot.cl[killi.annot.cl$Type %in% "CDS",]

##### Merge Info from killi annot and BLAST results
killi.mmu.rev <- merge(killi.annot.cl, blast.rev.res, by.x = "Nfur_Accession", by.y = "sseqid")

# select columns and reorder
colnames(killi.mmu.rev)
my.col.order <- c("Nfur_Accession", "Nfur_Description", "Nfur_Symbol", "Mmu_Symbol", "Mmu_ENSMUSP", "Mmu_ENSMUST", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore") 
killi.mmu.rev.cl <- killi.mmu.rev[,my.col.order]

write.table(killi.mmu.rev.cl, file = paste0(Sys.Date(),"_Mouse_Best_BLAST_hit_to_Killifish_Annotated_1e-3_withStats.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

##### Loop over table to retain only best hit
dim(killi.mmu.rev.cl) # [1] 63285    16

length(unique(killi.mmu.rev.cl$Nfur_Accession)) # 17969

killi.mmu.orth.rev <- unique(killi.mmu.rev.cl[,c("Nfur_Accession", "Nfur_Description", "Nfur_Symbol", "Mmu_Symbol") ])

write.table(killi.mmu.orth.rev, file = paste0(Sys.Date(),"_Mouse_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-3_Minimal_HOMOLOGY_TABLE_REV.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
