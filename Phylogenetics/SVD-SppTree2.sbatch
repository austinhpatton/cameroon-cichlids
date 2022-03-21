#!/bin/bash
#SBATCH --output=./err-out/SVDQ-SppTree2-NoOut_%J.out
#SBATCH --time=72:00:00
#SBATCH --account=fc_fishes
#SBATCH -p savio3

# infer a species tree using concatenated, linkage thinned SNPs using SVD quartets.
/global/scratch/users/austinhpatton/software/paup4a168_ubuntu64 -n /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/Phylogenies/BM-Ckot-WG-SNPs-SppTree2.min1.paup-block.txt
