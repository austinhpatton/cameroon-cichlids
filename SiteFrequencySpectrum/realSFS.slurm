#!/bin/bash
#SBATCH --output=./err-out/realSFS_%A_%a.out
#SBATCH --error=./err-out/realSFS_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --array=1-8
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=fc_fishes
#SBATCH -p savio3_bigmem

realSFS=/global/scratch/users/austinhpatton/software/angsd/misc/realSFS

id=$SLURM_ARRAY_TASK_ID

bamlist=$(ls bamlists/ | sed -n ${id}p)
pop=$(echo $bamlist | sed 's/-bams.txt//g')
outDir=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/angsd/saf/1d-saf-out

# find the global sfs
$realSFS -fold 1 -P 32 $outDir/$pop.saf.idx > /global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/angsd/saf/1d-sfs/$pop.sfs
