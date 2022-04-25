#!/bin/bash
#SBATCH --output=momi-sfs_%a.out
#SBATCH --error=momi-sfs_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=co_fishes
#SBATCH -p savio3

python=/global/scratch/users/austinhpatton/.conda/envs/momi/bin/python

vcfDir=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/FilteredThinnedSNPS

cd /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/momi/data

popfile=/global/home/users/austinhpatton/cichlids/cameroon/scripts/Onil_UMD/momi/FinalAnalysis/momi-popfile-final.txt
$python -m momi.read_vcf --no_aa --verbose --bed $vcfDir/snpable-mask-only-chroms.bed $vcfDir/momi-final-ChromsOnly-SNPs.vcf.gz $popfile BM-momi-FinalSamps.snpAlleleCounts.gz

$python -m momi.extract_sfs BM-momi-FinalSamps.sfs.gz 100 BM-momi-FinalSamps.snpAlleleCounts.gz
