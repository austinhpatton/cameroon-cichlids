#!/bin/bash

# This script is used by gnu parallel to build a BQSR model for each sample, to later be applied using the ApplyBQSR function in GATK. 

module load java

echo "Generating recalibration model for $1" 
/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx6G" BaseRecalibrator \
        -R /global/scratch/austinhpatton/cichlids/Oreochromis-Reference/Oreochromis.fna \
        -I /global/scratch/austinhpatton/cichlids/cameroon/preProcessing/mappedBams/$1_markdup.bam \
	--known-sites /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/Genotyped/bqsr1-BM-cichlids-filteredSNPS.vcf.gz \
	--known-sites /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/Genotyped/bqsr1-BM-cichlids_filtered_indels.vcf.gz \
        -O /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/RecalTables/$1_base-recal.table >& \
	/global/scratch/austinhpatton/cichlids/cameroon/gatk-output/RecalTables/$1-base-recal-table.out
echo "done"
