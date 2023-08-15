setwd('/Users/berenice/Dropbox/Manuscripts_and_Publications/2023/2023_Killifish_single_cell_Atlas_Bryan/TEST/Steroid_motif')
options(stringsAsFactors = F)

# 2023-06-08

#### Androgen responsive genes
are.homer <- read.csv('HOMER_ARE_annotations_to_genes_killi_Genome2015.xls', sep = "\t", header = T)
colnames(are.homer)[1] <- "PeakID"

### peaks within 2kb of TSS
hist(are.homer$Distance.to.TSS)

are.homer.2kbprom <- are.homer[which(abs(are.homer$Distance.to.TSS) < 2000),] # 5600
are.homer.1kbprom <- are.homer[which(abs(are.homer$Distance.to.TSS) < 1000),] # 3179
are.homer.500bprom <- are.homer[which(abs(are.homer$Distance.to.TSS) < 500),] # 1612

length(unique(are.homer.2kbprom$Gene.Name )) # 4795
length(unique(are.homer.1kbprom$Gene.Name )) # 2888
length(unique(are.homer.500bprom$Gene.Name)) # 1522

############################################################
#### Estrogen responsive genes
ere.homer <- read.csv('HOMER_ERE_annotations_to_genes_killi_Genome2015.xls', sep = "\t", header = T)
colnames(ere.homer)[1] <- "PeakID"

### peaks within 2kb of TSS
hist(ere.homer$Distance.to.TSS)

ere.homer.2kbprom <- ere.homer[which(abs(ere.homer$Distance.to.TSS) < 2000),] # 5148
ere.homer.1kbprom <- ere.homer[which(abs(ere.homer$Distance.to.TSS) < 1000),] # 2892
ere.homer.500bprom <- ere.homer[which(abs(ere.homer$Distance.to.TSS) < 500),] # 1516


length(unique(ere.homer.2kbprom$Gene.Name )) # 4536
length(unique(ere.homer.1kbprom$Gene.Name )) # 2679
length(unique(ere.homer.500bprom$Gene.Name)) # 1456

############
gmt.df <- data.frame(matrix(0,6,3))
gmt.df[,1] <- c("ARE_2kb"   ,
                "ARE_1kb"   ,
                "ARE_500bp" ,
                "ERE_2kb"   ,
                "ERE_1kb"   ,
                "ERE_500bp" )
gmt.df[,2] <- c("HOMER_ARE_2kb"   ,
                "HOMER_ARE_1kb"   ,
                "HOMER_ARE_500bp" ,
                "HOMER_ERE_2kb"   ,
                "HOMER_ERE_1kb"   ,
                "HOMER_ERE_500bp" )
gmt.df[,3] <- c(paste(unique(are.homer.2kbprom$Gene.Name  ), collapse = "\t"),
                paste(unique(are.homer.1kbprom$Gene.Name  ), collapse = "\t"),
                paste(unique(are.homer.500bprom$Gene.Name ), collapse = "\t"),
                paste(unique(ere.homer.2kbprom$Gene.Name  ), collapse = "\t"),
                paste(unique(ere.homer.1kbprom$Gene.Name  ), collapse = "\t"),
                paste(unique(ere.homer.500bprom$Gene.Name ), collapse = "\t"))


write.table(gmt.df, file = "2023-06-08_HOMER_2015_Killi_Genome_ERE_ARE_sites_Promoters.gmt", sep = "\t", quote = F, col.names = F, row.names = F)

            
            
            
            
            
            
