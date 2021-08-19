#!/bin/bash
module load java

echo "Running scaffold $1 for $2" 
/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx6G" HaplotypeCaller \
	-R /global/scratch/austinhpatton/cichlids/Oreochromis-Reference/Oreochromis.fna \
	-I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/BQSR-Bams/$2_BQSR.bam \
	--native-pair-hmm-threads $3 \
	-L $1 \
	-ERC GVCF \
	-O /global/scratch/austinhpatton/cichlids/cameroon/Post-BQSR/PerScaff/$2/$2-$1_raw_variants.g.vcf.gz >& \
	/global/scratch/austinhpatton/cichlids/cameroon/gatk-output/post-bqsr/$2/$2-$1-call-raw-variants.out
echo "done"
