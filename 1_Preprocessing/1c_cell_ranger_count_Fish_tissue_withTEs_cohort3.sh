~/Softwares/cellranger-6.1.2/bin/cellranger count --id Female_Kidney_3_FishTEDB_NR \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-01-24_killifish_tissues_cohort3/FASTQ \
                 --sample Female_Kidney_3_SIGAE5  \
                 --expect-cells 5000 \
                 --localmem 128 \
                 --localcores 16

~/Softwares/cellranger-6.1.2/bin/cellranger count --id Female_Liver_3_FishTEDB_NR \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-01-24_killifish_tissues_cohort3/FASTQ \
                 --sample Female_Liver_3_SIGAE4  \
                 --expect-cells 5000 \
                 --localmem 128 \
                 --localcores 16
                 
~/Softwares/cellranger-6.1.2/bin/cellranger count --id Female_Spleen_3_FishTEDB_NR \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-01-24_killifish_tissues_cohort3/FASTQ \
                 --sample Female_Spleen_3_SIGAE6  \
                 --expect-cells 5000 \
                 --localmem 128 \
                 --localcores 16

~/Softwares/cellranger-6.1.2/bin/cellranger count --id Male_Kidney_3_FishTEDB_NR \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-01-24_killifish_tissues_cohort3/FASTQ \
                 --sample Male_Kidney_3_SIGAE2  \
                 --expect-cells 5000 \
                 --localmem 128 \
                 --localcores 16

~/Softwares/cellranger-6.1.2/bin/cellranger count --id Male_Liver_3_FishTEDB_NR \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-01-24_killifish_tissues_cohort3/FASTQ \
                 --sample Male_Liver_3_SIGAE1  \
                 --expect-cells 5000 \
                 --localmem 128 \
                 --localcores 16

~/Softwares/cellranger-6.1.2/bin/cellranger count --id Male_Spleen_3_FishTEDB_NR \
                 --transcriptome ~/Softwares/cellranger-6.1.2/Reference/GCF_001465895.1_FishTEDB_NR \
                 --fastqs /mnt/data2/2022-01-24_killifish_tissues_cohort3/FASTQ \
                 --sample Male_Spleen_3_SIGAE3  \
                 --expect-cells 5000 \
                 --localmem 128 \
                 --localcores 16
