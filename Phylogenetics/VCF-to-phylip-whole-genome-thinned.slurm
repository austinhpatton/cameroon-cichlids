#!/bin/bash
#SBATCH --output=vcf2phylip_%J.out
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio2
#SBATCH --account=fc_fishes

module load python/3.7

cd /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/Phylogenies/
VCF=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/FilteredThinnedSNPS/BM-Ckott-LinkageThinned-FINAL-SNPs.vcf.gz

# Output a nexus file containing all sites across the genome. This may break the function if it's just too large. 
python /global/scratch/users/austinhpatton/software/vcf2phylip/vcf2phylip.py -i $VCF -n -f -p -m1 --output-folder /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/Phylogenies/ --output-prefix BM-Ckot-WG-SNPs-thinned

# dashes in the name confuses svdquartets
awk 'BEGIN {FS=OFS="|"}{gsub("-","_",$1)}1' /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/Phylogenies/BM-Ckot-WG-SNPs-thinned.min1.nexus > tmp && mv tmp /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/Phylogenies/BM-Ckot-WG-SNPs-thinned.min1.nexus
