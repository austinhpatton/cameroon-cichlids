#!/bin/bash
module load java

# This script is to be submitted via a gnu-parallel submission (included in this repository). 
# This conducted the first round of genotyping to build our base score quality recalibration model with.

# $1 is the scaffold used as input, $2 is the sample, $3 is the number of threads used, as determined 
# in the associated slurm script. 

echo "Running scaffold $1 for $2"
/global/scratch/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx6G" HaplotypeCaller \
        -R /global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fa \
        -I /global/scratch/users/austinhpatton/cichlids/cameroon/preProcessing/mappedBams_OnilUMD/$1_markdup.bam \
        --native-pair-hmm-threads $2 \
        -L /global/home/users/austinhpatton/cichlids/cameroon/scripts/Onil_UMD/GATK/PreBQSR/ScaffNames.list \
        -ERC GVCF \
        -O /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/$1/$1-ShortScaffs_raw_variants.g.vcf.gz >& \
        /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/gatk-output/pre-bqsr/$1/$1-ShortScaffs-call-raw-variants.out
echo "done"
