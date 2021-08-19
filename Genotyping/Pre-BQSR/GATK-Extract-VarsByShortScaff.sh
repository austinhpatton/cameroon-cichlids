#!/bin/bash

# This is the shell script invoked by the parallel slurm script for 
# calculating genotype likelihoods from all remaining short scaffolds.

module load java

/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx6G" HaplotypeCaller \
        -R /global/scratch/austinhpatton/cichlids/Oreochromis-Reference/Oreochromis.fna \
        -I /global/scratch/austinhpatton/cichlids/cameroon/preProcessing/mappedBams/$1_markdup.bam \
	--native-pair-hmm-threads $2 \
        -L ScaffNames.list \
        -ERC GVCF \
        -O /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-ShortScaffs_raw_variants.g.vcf.gz >& \
	/global/scratch/austinhpatton/cichlids/cameroon/gatk-output/bqsr-prep/bqsr1/$1/$1-ShortScaffs-call-raw-variants.out
