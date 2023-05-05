WORKDIR=/scratch2/alanxu/dietary_liver/align
cd $WORKDIR
/home1/alanxu/subread-2.0.4-Linux-x86_64/bin/featureCounts -a /scratch2/alanxu/2019_killi_rna_raw/deconv/genes.gtf -F GTF -O -T 20 -p --countReadPairs -o /scratch2/alanxu/dietary_liver/count/liver.txt \
*.sortedByCoord.out.bam