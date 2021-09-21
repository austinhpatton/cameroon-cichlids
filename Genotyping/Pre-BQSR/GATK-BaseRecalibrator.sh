#!/bin/bash

# This script is used by gnu parallel to build a BQSR model for each sample, to later be applied using the ApplyBQSR function in GATK. 

module load java

echo "Generating recalibration model for $1" 
/global/scratch/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx6G" BaseRecalibrator \
        -R /global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fa \
        -I /global/scratch/users/austinhpatton/cichlids/cameroon/preProcessing/mappedBams_OnilUMD/$1_markdup.bam \
	--known-sites /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/HardFilteredVariants/FINAL-BM-cichlids-filteredSNPS.vcf.gz \
	--known-sites /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/HardFilteredVariants/FINAL-BM-cichlids_filtered_indels.vcf.gz \
        -O /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/RecalTables/$1_base-recal.table >& \
	/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/RecalTables/$1-base-recal-table.out
echo "done"
