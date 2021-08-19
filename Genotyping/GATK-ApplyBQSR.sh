#!/bin/bash

# This script applies the fitted BQSR model to each sample, in parallel. This is invoked by "GATK-ApplyBQSR-Parallel.slurm"

module load java

echo "Applying BQSR to $1" 
/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx6G" ApplyBQSR \
  -R /global/scratch/austinhpatton/cichlids/Oreochromis-Reference/Oreochromis.fna \
  -I /global/scratch/austinhpatton/cichlids/cameroon/preProcessing/mappedBams/$1_markdup.bam \
  --bqsr-recal-file /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/RecalTables/$1_base-recal.table \
  -O /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/BQSR-Bams/$1_BQSR.bam >& /global/scratch/austinhpatton/cichlids/cameroon/gatk-output/ApplyBQSR/$1-apply-BQSR.out
echo "done"
