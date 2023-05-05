setwd('/Volumes/BB_Home_HQ2/SIngle_Cell_analysis/2022-02-01_Killifish_scRNAseq_ATLAS_analyses/Zebrafish_alignment')
options(stringsAsFactors = FALSE)

#########################################################################################################################
# read in result table and add in column names
blast.rev.res <- read.table('2022-02-22_Zebrafish_Nfur2015_BestHits_1e-5_REVERSE.txt', sep = "\t", header = F)
colnames(blast.rev.res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast.rev.res$DanRer_Symbol  <- unlist(lapply(strsplit(blast.rev.res$qseqid,"|", fixed = TRUE),'[',1))
blast.rev.res$DanRer_ENSDARP <- unlist(lapply(strsplit(blast.rev.res$qseqid,"|", fixed = TRUE),'[',2))
blast.rev.res$DanRer_ENSDART <- unlist(lapply(strsplit(blast.rev.res$qseqid,"|", fixed = TRUE),'[',3))

# read in killifish gene annotation table
killi.annot    <- read.csv('GCF_001465895.1_Nfu_20140520_feature_table.txt', header = F, sep = "\t")
killi.annot.cl <- killi.annot[,c(1,11,14, 15)]
colnames(killi.annot.cl) <- c("Type", "Nfur_Accession", "Nfur_Description", "Nfur_Symbol")

# select on protein only
killi.annot.cl <- killi.annot.cl[killi.annot.cl$Type %in% "CDS",]


##### Merge Info from killi annot and BLAST results
killi.zebra.rev <- merge(killi.annot.cl, blast.rev.res, by.x = "Nfur_Accession", by.y = "sseqid")

# select columns and reorder
colnames(killi.zebra.rev)
my.col.order <- c("Nfur_Accession", "Nfur_Description", "Nfur_Symbol", "DanRer_Symbol", "DanRer_ENSDARP", "DanRer_ENSDART", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore") 
killi.zebra.rev.cl <- killi.zebra.rev[,my.col.order]

write.table(killi.zebra.rev.cl, file = paste0(Sys.Date(),"_Zebrafish_Best_BLAST_hit_to_Killifish_Annotated_1e-5_withStats.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

##### Loop over table to retain only best hit
dim(killi.zebra.rev.cl) # [1] 38942    16

length(unique(killi.zebra.rev.cl$Nfur_Accession)) # 19503

killi.zebra.orth.rev <- unique(killi.zebra.rev.cl[,c("Nfur_Accession", "Nfur_Description", "Nfur_Symbol", "DanRer_Symbol", "DanRer_ENSDARP", "DanRer_ENSDART") ])

write.table(killi.zebra.orth.rev, file = paste0(Sys.Date(),"_Zebrafish_Best_BLAST_hit_to_Killifish_Annotated_hit_1e-5_Minimal_HOMOLOGY_TABLE_REV.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
