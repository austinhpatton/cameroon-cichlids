#!/bin/bash
#SBATCH --output=genotype-BSQR1_%a.out
#SBATCH --error=genotype-BSQR1_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio3
#SBATCH --account=co_fishes

# This script combines the per-sample g.vcfs, and subsequently conducts joint genotyping on the cohort. 
# Then, SNPs and Indels are pulled out, and subsequently hard filtered. These will be the subject of 
# Base Score Quality Recalibration (BQSR). 

module load python
module load samtools
module load bcftools
module load bwa 
module load java

REF='/global/scratch/austinhpatton/cichlids/Oreochromis-Reference/Oreochromis.fna'
INTERIM='/global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/Genotyped'
bamDir='/global/scratch/austinhpatton/cichlids/cameroon/preProcessing/mappedBams/'
gatkOut='/global/scratch/austinhpatton/cichlids/cameroon/gatk-output/bqsr-prep/'

gatk="/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk"

PerSampDir="/global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/PerSamp"

# Combine GVCF files produced by HaplotypeCaller
$gatk --java-options "-Xmx192g" CombineGVCFs \
	-R $REF -O $INTERIM/bqsr1-cichlids-combined.g.vcf.gz \
	-V $PerSampDir/CAR152_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/CAR389_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/DIK101_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/DIK102_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/DIK103_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/DIK104_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/DIK201_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/DIK208_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/DIK240_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/DIK295_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/DIK352_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/DIK353_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/EIS072_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/EIS084_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/EIS085_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/EIS100_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/EIS101_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/EIS102_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/EIS103_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/EIS104_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/EIS105_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/EIS227_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/KON011_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/KON028_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/Kdik-HighCov-AllSites.g.vcf.gz \
	-V $PerSampDir/Kdik248_L005-AllSites.g.vcf.gz \
	-V $PerSampDir/Kdik295_L006-AllSites.g.vcf.gz \
	-V $PerSampDir/Keis-HighCov-AllSites.g.vcf.gz \
	-V $PerSampDir/Keis332_L003-AllSites.g.vcf.gz \
	-V $PerSampDir/MA10_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MA11_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MA14_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MA3_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MA5_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MA6_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MA7_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MA8_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MA9_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MAR038_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MAR092_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MAR101_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MAR104_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MAR105_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MAR106_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MAR107_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MAR108_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/ME4_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MON012_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON013_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON015_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON016_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON021_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON100_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON103_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON104_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON105_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON106_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON107_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON108_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MON109_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/MU3_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MU5_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MU6_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MU7_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/MYA252_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/Mmya-HighCov-AllSites.g.vcf.gz \
	-V $PerSampDir/Mmya210_L006-AllSites.g.vcf.gz \
	-V $PerSampDir/Mmya252_L005-AllSites.g.vcf.gz \
	-V $PerSampDir/PIN02_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PIN040_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PIN050_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PIN071_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PIN075_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PIN083_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PIN100_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PIN101_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PIN102_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PIN103_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PUN100_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/PUN101_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/PUN102_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/PUN104_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/PUN105_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/PUN107_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/PUN108_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/PUN109_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/PUN110_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/Pmac-HighCov-AllSites.g.vcf.gz \
	-V $PerSampDir/Pmac226_L006-AllSites.g.vcf.gz \
	-V $PerSampDir/Pmac276_L006-AllSites.g.vcf.gz \
	-V $PerSampDir/STE102_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/STE106_L004-AllSites.g.vcf.gz \
	-V $PerSampDir/STE107_L001-AllSites.g.vcf.gz \
	-V $PerSampDir/Scar-HighCov-AllSites.g.vcf.gz \
	-V $PerSampDir/Scar152_L003-AllSites.g.vcf.gz \
	-V $PerSampDir/Scar226_L003-AllSites.g.vcf.gz \
	-V $PerSampDir/Scar389_L003-AllSites.g.vcf.gz \
	-V $PerSampDir/Sgal-CRM-HighCov-AllSites.g.vcf.gz \
	-V $PerSampDir/Sgal-Meme-HighCov-AllSites.g.vcf.gz \
	-V $PerSampDir/Sgal-Mongo-HighCov-AllSites.g.vcf.gz \
	-V $PerSampDir/SgalME1_L005-AllSites.g.vcf.gz \
	-V $PerSampDir/SgalME2_L006-AllSites.g.vcf.gz \
	-V $PerSampDir/SgalMU1_L006-AllSites.g.vcf.gz \
	-V $PerSampDir/SgalMU2_L005-AllSites.g.vcf.gz \
	-V $PerSampDir/Skna001_L002-AllSites.g.vcf.gz \
	-V $PerSampDir/Skna053_L002-AllSites.g.vcf.gz \
	-V $PerSampDir/Slam225_L002-AllSites.g.vcf.gz \
	-V $PerSampDir/Slin220_L003-AllSites.g.vcf.gz \
	-V $PerSampDir/Sloh480_L002-AllSites.g.vcf.gz \
	-V $PerSampDir/Smar-HighCov-AllSites.g.vcf.gz \
	-V $PerSampDir/Smar231_L005-AllSites.g.vcf.gz \
	-V $PerSampDir/SmarNOL_L006-AllSites.g.vcf.gz \
	-V $PerSampDir/Smon-HighCov-AllSites.g.vcf.gz \
	-V $PerSampDir/Smon001_L005-AllSites.g.vcf.gz \
	-V $PerSampDir/Spin-10X_L007-AllSites.g.vcf.gz \
	-V $PerSampDir/Spin332_L002-AllSites.g.vcf.gz \
	-V $PerSampDir/SpinNOL_L005-AllSites.g.vcf.gz \
	-V $PerSampDir/TguiMA1_L005-AllSites.g.vcf.gz \
	-V $PerSampDir/TguiMA2_L006-AllSites.g.vcf.gz \
	-V $PerSampDir/TguiMA4_L005-AllSites.g.vcf.gz 

# Perform joint genotyping on combined GVCF
echo "Performing joint genotyping"
$gatk --java-options "-Xmx192g" GenotypeGVCFs -R $REF -V $INTERIM/bqsr1-BM-cichlids-combined.g.vcf.gz -O $INTERIM/bqsr1-joint-gt-BM-cichlids.vcf.gz 

echo "Pulling out SNPs"
$gatk SelectVariants -V $INTERIM/bqsr1-joint-gt-BM-cichlids.vcf.gz -select-type SNP -O $INTERIM/BM-cichlids_bqsr1_snps_init.vcf.gz >& $INTERIM/BM-cichlids-bqsr1-snps.out

echo "Pulling out INDELs"
$gatk SelectVariants -V $INTERIM/bqsr1-joint-gt-BM-cichlids.vcf.gz -select-type INDEL -O $INTERIM/BM-cichlids_bqsr1_indels_init.vcf.gz >& $INTERIM/BM-cichlids-bqsr1-indels.out

# Hard filter the SNPs
echo "Hard filtering SNPs"
$gatk VariantFiltration -R $REF -V $INTERIM/BM-cichlids_bqsr1_snps_init.vcf.gz -O $INTERIM/bqsr1-BM-cichlids-filteredSNPS.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" >& $INTERIM/BM-cichlids-bqsr1-filter-raw-snps.out

# Hard filter the indels
echo "Hard filtering INDELs"
$gatk VariantFiltration -R $REF -V $INTERIM/BM-cichlids_bqsr1_indels_init.vcf.gz -O $INTERIM/bqsr1-BM-cichlids_filtered_indels.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" \
	-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" >& $INTERIM/BM-cichlids-bqsr1-filter-raw-indels.out
