#!/bin/bash
#SBATCH --output=dbimport_%A_%a.out
#SBATCH --error=dbimport_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --array=1-23
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=fc_fishes
#SBATCH -p savio3

id=$SLURM_ARRAY_TASK_ID
scaff=$(sed "${id}q;d" regions.txt)
region=$(sed "${id}q;d" regions.txt)

if [ $id == 23 ]
then
	region=ScaffNames.list
	scaff=ShortScaffs
else
	region=$(echo $scaff)	
fi

module load java

/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx5G" GenomicsDBImport -L ${region} \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/CAR152_L004/CAR152_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/CAR389_L001/CAR389_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/DIK101_L004/DIK101_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/DIK102_L001/DIK102_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/DIK103_L001/DIK103_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/DIK104_L001/DIK104_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/DIK201_L004/DIK201_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/DIK208_L001/DIK208_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/DIK240_L001/DIK240_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/DIK295_L001/DIK295_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/DIK352_L001/DIK352_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/DIK353_L001/DIK353_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/EIS072_L001/EIS072_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/EIS084_L004/EIS084_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/EIS085_L001/EIS085_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/EIS100_L001/EIS100_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/EIS101_L001/EIS101_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/EIS102_L001/EIS102_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/EIS103_L001/EIS103_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/EIS104_L001/EIS104_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/EIS105_L001/EIS105_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/EIS227_L001/EIS227_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/KON011_L001/KON011_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/KON028_L001/KON028_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Kdik-HighCov/Kdik-HighCov-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Kdik248_L005/Kdik248_L005-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Kdik295_L006/Kdik295_L006-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Keis-HighCov/Keis-HighCov-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Keis332_L003/Keis332_L003-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MA10_L001/MA10_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MA11_L001/MA11_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MA14_L001/MA14_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MA3_L001/MA3_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MA5_L001/MA5_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MA6_L001/MA6_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MA7_L001/MA7_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MA8_L001/MA8_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MA9_L001/MA9_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MAR038_L004/MAR038_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MAR092_L004/MAR092_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MAR101_L004/MAR101_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MAR104_L001/MAR104_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MAR105_L001/MAR105_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MAR106_L001/MAR106_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MAR107_L001/MAR107_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MAR108_L001/MAR108_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/ME4_L001/ME4_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON012_L004/MON012_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON013_L004/MON013_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON015_L004/MON015_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON016_L004/MON016_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON021_L004/MON021_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON100_L004/MON100_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON103_L004/MON103_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON104_L004/MON104_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON105_L004/MON105_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON106_L004/MON106_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON107_L004/MON107_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON108_L004/MON108_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MON109_L004/MON109_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MU3_L001/MU3_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MU5_L001/MU5_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MU6_L001/MU6_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MU7_L001/MU7_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/MYA252_L001/MYA252_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Mmya-HighCov/Mmya-HighCov-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Mmya210_L006/Mmya210_L006-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Mmya252_L005/Mmya252_L005-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PIN02_L004/PIN02_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PIN040_L004/PIN040_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PIN050_L004/PIN050_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PIN071_L004/PIN071_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PIN075_L004/PIN075_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PIN083_L004/PIN083_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PIN100_L004/PIN100_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PIN101_L004/PIN101_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PIN102_L004/PIN102_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PIN103_L004/PIN103_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PUN100_L001/PUN100_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PUN101_L001/PUN101_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PUN102_L001/PUN102_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PUN104_L004/PUN104_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PUN105_L001/PUN105_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PUN107_L001/PUN107_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PUN108_L001/PUN108_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PUN109_L001/PUN109_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/PUN110_L001/PUN110_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Pmac-HighCov/Pmac-HighCov-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Pmac226_L006/Pmac226_L006-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Pmac276_L006/Pmac276_L006-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/STE102_L001/STE102_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/STE106_L004/STE106_L004-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/STE107_L001/STE107_L001-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Scar-HighCov/Scar-HighCov-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Scar152_L003/Scar152_L003-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Scar226_L003/Scar226_L003-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Scar389_L003/Scar389_L003-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Sgal-CRM-HighCov/Sgal-CRM-HighCov-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Sgal-Meme-HighCov/Sgal-Meme-HighCov-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Sgal-Mongo-HighCov/Sgal-Mongo-HighCov-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/SgalME1_L005/SgalME1_L005-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/SgalME2_L006/SgalME2_L006-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/SgalMU1_L006/SgalMU1_L006-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/SgalMU2_L005/SgalMU2_L005-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Skna001_L002/Skna001_L002-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Skna053_L002/Skna053_L002-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Slam225_L002/Slam225_L002-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Slin220_L003/Slin220_L003-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Sloh480_L002/Sloh480_L002-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Smar-HighCov/Smar-HighCov-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Smar231_L005/Smar231_L005-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/SmarNOL_L006/SmarNOL_L006-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Smon-HighCov/Smon-HighCov-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Smon001_L005/Smon001_L005-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Spin-10X_L007/Spin-10X_L007-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Spin332_L002/Spin332_L002-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/SpinNOL_L005/SpinNOL_L005-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/TguiMA1_L005/TguiMA1_L005-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/TguiMA2_L006/TguiMA2_L006-${scaff}_raw_variants.g.vcf.gz \
	-V /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/TguiMA4_L005/TguiMA4_L005-${scaff}_raw_variants.g.vcf.gz \
	--genomicsdb-workspace-path /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Joint/${scaff} \
	--tmp-dir /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Joint/TMP
