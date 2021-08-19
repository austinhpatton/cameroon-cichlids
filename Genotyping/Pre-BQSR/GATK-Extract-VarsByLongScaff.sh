#!/bin/bash
module load java

# This script is to be submitted via a gnu-parallel submission (included in this repository). 
# This conducted the first round of genotyping to build our base score quality recalibration model with.

# $1 is the scaffold used as input, $2 is the sample, $3 is the number of threads used, as determined 
# in the associated slurm script. 

echo "Running scaffold $1 for $2"
/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx6G" HaplotypeCaller \
        -R /global/scratch/austinhpatton/cichlids/Oreochromis-Reference/Oreochromis.fna \
        -I /global/scratch/austinhpatton/cichlids/cameroon/preProcessing/mappedBams/$2_markdup.bam \
        --native-pair-hmm-threads $3 \
        -L $1 \
        -ERC GVCF \
        -O /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$2/$2-$1_raw_variants.g.vcf.gz >& \
        /global/scratch/austinhpatton/cichlids/cameroon/gatk-output/bqsr-prep/bqsr1/$2/$2-$1-call-raw-variants.out
echo "done"
