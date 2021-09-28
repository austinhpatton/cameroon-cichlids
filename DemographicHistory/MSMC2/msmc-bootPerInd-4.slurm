#!/bin/bash
#SBATCH --output=./err-out/bs.msmc_%A_%a.out
#SBATCH --error=./err-out/bs.msmc_%A_%a.err
#SBATCH --time=72:00:00
#SBATCH --array=1-11
#SBATCH --ntasks-per-node=1
#SBATCH --account=co_rosalind
###SBATCH --cpus-per-task=15
###SBATCH --mem-per-cpu=2G
#SBATCH -p savio3

module load msmc2
module load python/3.7

MSMCTOOLS="/global/scratch/users/austinhpatton/software/msmc-tools"
msmcIn=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/MSMC/Input
bootDir=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/MSMC/Bootstrap

id=$SLURM_ARRAY_TASK_ID
samp=$(sed "${id}q;d" msmc-samps.txt)

# Set the input for bootstrapping
bootIn=$(for scaff in $(cat ChromNames.txt); do find $msmcIn/ -maxdepth 1 -name "${samp}-${scaff}-msmc-multihetsep.txt"; done)

# And set output
bootOut=$bootDir/${samp}.bootstrap

echo "generating bootstraps for ${POP_OR_IND}.${RUN_NAME}"
$MSMCTOOLS/multihetsep_bootstrap.py -n 50 -s 2272727 --chunks_per_chromosome 20 --nr_chromosomes 22 $bootOut $bootIn

# Now, run MSMC, looping through bootstraps
for i in {1..50}
do
	# Now, run MSMC
	# Set a meaningful output name
	msmcOut="${samp}_bootstrap_${i}"
	echo "Beginning bootstrapping procedure for ${samp}, iteration $i"
	~/software/msmc2-2.1.3/build/release/msmc2 -t 30 -p 1*2+18*1+1*2+1*3 -i 50 -I 0,1 -o $bootDir/${samp}.bootstrap_${i}/$msmcOut $bootDir/${samp}*_${i}/*
done

