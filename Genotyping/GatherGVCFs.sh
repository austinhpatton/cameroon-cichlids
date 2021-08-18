#!/bin/bash

module load java

/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx5G" GatherVcfs \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022199.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022200.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022201.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022202.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022203.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022204.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022205.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022206.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022207.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022208.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022209.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022210.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022211.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022212.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022213.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022214.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022215.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022216.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022217.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022218.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022219.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-NC_022220.1_raw_variants.g.vcf.gz \
        -I /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerScaff/$1/$1-ShortScaffs_raw_variants.g.vcf.gz \
        --CREATE_INDEX false \
        -O /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerSamp/$1-AllSites.g.vcf

~/software/htslib-1.11/bgzip /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerSamp/$1-AllSites.g.vcf
~/software/htslib-1.11/tabix -p vcf /global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerSamp/$1-AllSites.g.vcf.gz
