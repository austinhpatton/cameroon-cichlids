#!/bin/bash
#SBATCH --output=iqtree_%J.out
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=fc_fishes
#SBATCH -p savio3

/global/scratch/users/austinhpatton/software/iqtree-1.6.12-Linux/bin/iqtree -o Ckot383_L002 -o Ckot499_L002 -st DNA -m GTR+ASC -nt 32 -alrt 1000 -bb 1000 -s /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/Phylogenies/BM-Ckot-WG-SNPs-thinned.min1.fasta.varsites.phy
