#!/bin/bash
#SBATCH --output=filt-genotypes_%a.out
#SBATCH --error=filt-genotypes_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio3
#SBATCH --account=fc_fishes

module load java

gatk=/global/scratch/users/austinhpatton/software/gatk-4.1.8.1/gatk
REF=/global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fa
BASE=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/Genotyped
FILT=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/HardFilteredVariants

# Now we will select variants, and conduct hard filtering with GATK. We then will use these as the reference "high quality" panel of variants for BQSR.

echo "Pulling out SNPs"
$gatk SelectVariants -V $BASE/all-BM-cichlids-joint-gt-post-bqsr-RAW.vcf.gz -select-type SNP -O $BASE/BM-cichlids_post-bqsr_snps_init.vcf.gz >& $BASE/BM-cichlids-post-bqsr-snps.out

echo "Pulling out INDELs"
$gatk SelectVariants -V $BASE/all-BM-cichlids-joint-gt-post-bqsr-RAW.vcf.gz -select-type INDEL -O $BASE/BM-cichlids_post-bqsr_indels_init.vcf.gz >& $BASE/BM-cichlids-post-bqsr-indels.out

echo "Hard filtering SNPs"
# Hard filter the SNPs
$gatk VariantFiltration -R $REF -V $BASE/BM-cichlids_post-bqsr_snps_init.vcf.gz -O $FILT/FINAL-BM-cichlids-filteredSNPS.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" >& $BASE/BM-cichlids-post-bqsr-filter-raw-snps.out

echo "Hard filtering INDELs"
# Hard filter the indels
$gatk VariantFiltration -R $REF -V $BASE/BM-cichlids_post-bqsr_indels_init.vcf.gz -O $FILT/FINAL-BM-cichlids_filtered_indels.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" >& $BASE/BM-cichlids-post-bqsr-filter-raw-indels.out

