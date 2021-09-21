#!/bin/bash

# This script applies the fitted BQSR model to each sample, in parallel. This is invoked by "GATK-ApplyBQSR-Parallel.slurm"

module load java

echo "Applying BQSR to $1" 
/global/scratch/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx6G" ApplyBQSR \
  -R /global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fa \
  -I /global/scratch/users/austinhpatton/cichlids/cameroon/preProcessing/mappedBams_OnilUMD/$1_markdup.bam \
  --bqsr-recal-file /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/RecalTables/$1_base-recal.table \
  -O /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/RecalBams/$1_BQSR.bam >& /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/gatk-output/ApplyBQSR/$1-apply-BQSR.out
echo "done"
