# Author: Alan Xu
# This script uses star/2.7.0e to align polished reads

module purge
module load gcc/8.3.0
module load intel/18.0.4
module load star/2.7.0e

WORKDIR=/scratch2/alanxu/dietary_liver
cd $WORKDIR
for f in $(find . -name '*_1_trimmed.fastq')
do
    f1=$(basename "${f}")
    f2=$(basename "${f}" | sed 's/_1_trimmed.fastq/_2_trimmed.fastq/')
    outname="$(echo $f1 | sed 's/_1_trimmed.fastq//')"
    STAR --genomeDir /scratch2/alanxu/2019_killi_rna_raw/deconv/2015_genome_index --readFilesIn $f1 $f2 \
    --runThreadN 12 --outFilterMultimapNmax 200 --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --alignEndsProtrude 10 ConcordantPair --limitGenomeGenerateRAM 60000000000 --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix /scratch2/alanxu/dietary_liver/count/$outname
done