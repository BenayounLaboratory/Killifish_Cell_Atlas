# Author: Alan Xu
# Software Used: fastp version 0.23.2

WORKDIR=/scratch2/alanxu/dietary_liver
cd $WORKDIR
for f in $(find . -name '*_1.fastq')
do
    f1=$(basename "${f}")
    f2=$(basename "${f}" | sed 's/_1.fastq/_2.fastq/')
    o1="$(echo $f1 | sed 's/.fastq/_trimmed.fastq/')"
    o2="$(echo $f2 | sed 's/.fastq/_trimmed.fastq/')"
    ~/fastp --in1 $f1 --in2 $f2 --out1 $o1 --out2 $o2 --thread 20
done