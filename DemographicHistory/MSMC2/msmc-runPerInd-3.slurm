#!/bin/bash
#SBATCH --output=./err-out/run.msmc_%A_%a.out
#SBATCH --error=./err-out/run.msmc_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --array=1-11
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=fc_fishes
###SBATCH --mem=348G
###SBATCH --cpus-per-task=15
###SBATCH --mem-per-cpu=2G
#SBATCH -p savio3

module load msmc2

msmcIn=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/MSMC/Input
msmcOut=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/MSMC/Output
msmc2=/global/scratch/users/austinhpatton/software/msmc2-2.1.3/build/release/msmc2

id=$SLURM_ARRAY_TASK_ID
samp=$(sed "${id}q;d" msmc-samps.txt)

# This script will run an MSMC analysis per-individual via an array job

# First, list all per-scaffold multihetsep files per individual to a file
ls $msmcIn/${samp}*msmc-multihetsep.txt > $msmcIn/${samp}-AllInput

# And concatenate these into a single file to estimatejoint demography from all chromosomes
for scaff in $(cat $msmcIn/${samp}-AllInput)
do 
	cat $scaff >> $msmcIn/${samp}-AllInput-multihetsep.txt
done

# Now run MSMC2
$msmc2 -t 30 -p 1*2+22*1+1*2+1*2 -i 50 -o $msmcOut/$samp -I 0,1 $msmcIn/${samp}*msmc-multihetsep.txt
