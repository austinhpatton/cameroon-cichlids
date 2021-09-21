#!/bin/bash
#SBATCH --output=./err-out/JointGenotypeScaffs_%a.out
#SBATCH --error=./err-out/JointGenotypeScaffs_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --array=1-23
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=fc_fishes
#SBATCH -p savio3

module load python
module load samtools
module load bcftools
module load bwa 
module load java

REF='/global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fa'
FINAL='/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/PreBQSR'
bamDir='/global/scratch/users/austinhpatton/cichlids/cameroon/preProcessing/mappedBams_OnilUMD/'
gatkOut='/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/gatk-output/pre-bqsr'
OUT='/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/PreBQSR/Genotyped'

id=$SLURM_ARRAY_TASK_ID
scaff=$(sed "${id}q;d" DBs-To-Genotype.txt)

# Combine GVCF files produced by HaplotypeCaller
## Perform joint genotyping on per-scaff genomics DB - much more efficient this way. 
echo "Conducting Joint Genotyping for $scaff"
cd /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Joint/
$gatk --java-options "-Xmx192g" GenotypeGVCFs \
    -R $REF \
    -V gendb://${scaff} \
    -O $OUT/${scaff}-joint-gt-BM-cichlids-post-bqsr.vcf.gz
