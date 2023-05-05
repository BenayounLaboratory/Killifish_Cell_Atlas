#convert_cell_ranger_to_loom.sh

#this script uses velocyto0.17.17 to convert cell ranger out folders to loom files

gtf="tissue_atlas_gtf.gtf"

for f in $(find "." -name '*_FishTEDB_NR')

do
velocyto run10x -@ 10 --samtools-memory 1000 $f $gtf  
done
