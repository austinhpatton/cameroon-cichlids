#!/bin/bash
module load java

/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx6G" HaplotypeCaller \
	-R /global/scratch/austinhpatton/cichlids/Oreochromis-Reference/Oreochromis.fna \
	-I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/BQSR-Bams/$1_BQSR.bam \
	--native-pair-hmm-threads $2 \
	-L /global/home/users/austinhpatton/cichlids/cameroon/scripts/StreamlinedGATK/Post-BQSR/ScaffNames.list \
	-ERC GVCF \
	-O /global/scratch/austinhpatton/cichlids/cameroon/Post-BQSR/PerScaff/$1/$1-ShortScaffs_raw_variants.g.vcf.gz >& \
	/global/scratch/austinhpatton/cichlids/cameroon/gatk-output/post-bqsr/$1/$1-ShortScaffs-call-raw-variants.out
