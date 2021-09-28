#!/bin/bash
#SBATCH --output=./err-out/1d-saf_%A_%a.out
#SBATCH --error=./err-out/1d-saf_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --array=1-8
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=co_rosalind
#SBATCH -p savio3

angsd=/global/scratch/users/austinhpatton/software/angsd/angsd
realSFS=/global/scratch/users/austinhpatton/software/angsd/misc/realSFS

id=$SLURM_ARRAY_TASK_ID

bamlist=$(ls bamlists/ | sed -n ${id}p)
pop=$(echo $bamlist | sed 's/-bams.txt//g')

ref=/global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fa
outDir=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/angsd/saf/1d-saf-out

# Calculate the site allele frequency spectrum from genotype likelihoods
$angsd -b ./bamlists/$bamlist -anc $ref -out $outDir/$pop -P 32 -dosaf 1 -gl 1 -minMapQ 30 -minQ 20

# find the global sfs
$realSFS -fold 1 -P 32 $outDir/$pop.saf.idx
