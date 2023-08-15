scanMotifGenomeWide.pl are.motif 2021-11-30_GCF_001465895.1_MASKED_withFishTEDBcontigs.fa -bed > 2023-06-08_HOMER_ARE.sites.killi2015_genome.bed
scanMotifGenomeWide.pl ere.motif 2021-11-30_GCF_001465895.1_MASKED_withFishTEDBcontigs.fa -bed > 2023-06-08_HOMER_ERE.sites.killi2015_genome.bed

##    39426 2023-06-08_HOMER_ARE.sites.killi2015_genome.bed
##    37637 2023-06-08_HOMER_ERE.sites.killi2015_genome.bed

#### cat 2023-06-08_HOMER_ARE.sites.killi2015_genome.bed | perl -lane '@bla = split(/\t/,$_); @bla2 = split(/\s/,$bla[0]); print "$bla2[0]\t$bla[1]\t$bla[2]\t$bla[3]\t$bla[4]"' > 2023-06-08_HOMER_ARE.sites.killi2015_genome.CLEAN.bed
#### cat 2023-06-08_HOMER_ERE.sites.killi2015_genome.bed | perl -lane '@bla = split(/\t/,$_); @bla2 = split(/\s/,$bla[0]); print "$bla2[0]\t$bla[1]\t$bla[2]\t$bla[3]\t$bla[4]"' > 2023-06-08_HOMER_ERE.sites.killi2015_genome.CLEAN.bed

annotatePeaks.pl 2023-06-08_HOMER_ARE.sites.killi2015_genome.CLEAN.bed 2021-11-30_GCF_001465895.1_MASKED_withFishTEDBcontigs.fa -gtf GCF_001465895.1_Nfu_20140520_genomic_CLEAN_MT_exon.filtered.gtf > HOMER_ARE_annotations_to_genes_killi_Genome2015.xls
annotatePeaks.pl 2023-06-08_HOMER_ERE.sites.killi2015_genome.CLEAN.bed 2021-11-30_GCF_001465895.1_MASKED_withFishTEDBcontigs.fa -gtf GCF_001465895.1_Nfu_20140520_genomic_CLEAN_MT_exon.filtered.gtf > HOMER_ERE_annotations_to_genes_killi_Genome2015.xls