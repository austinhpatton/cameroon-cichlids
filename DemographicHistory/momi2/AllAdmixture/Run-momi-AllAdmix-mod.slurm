#!/bin/bash
#SBATCH --output=AllAdmix-fit-cont.out
#SBATCH --error=AllAdmix-fit-cont.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=fc_fishes
#SBATCH -p savio3_bigmem

python=/global/scratch/users/austinhpatton/.conda/envs/momi/bin/python
jupyter=/global/scratch/users/austinhpatton/.conda/envs/momi/bin/jupyter

conda init momi

# The timeout flag is so that the jupyter notebook doesn't prematurely stop since it's being run on the cluster. 
$jupyter nbconvert --ExecutePreprocessor.timeout=999999999999 --execute --to html Bmbo-momi-AllAdmix-InitFit.ipynb 
