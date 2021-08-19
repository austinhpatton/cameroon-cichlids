#!/bin/bash
#SBATCH --output=ExtractVars_set-00_%a.out
#SBATCH --error=ExtractVars_set-00_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --account=ac_fishes
#SBATCH -p savio3

# Load the default version of GNU parallel.
module load gnu-parallel
module load java 

REF='/global/scratch/austinhpatton/cichlids/Oreochromis-Reference/Oreochromis.fna'
FINAL='/global/scratch/austinhpatton/cichlids/cameroon/Post-BQSR/'
bamDir='/global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/bqsr1/BQSR-Bams/'
gatkOut='/global/scratch/austinhpatton/cichlids/cameroon/gatk-output/post-bqsr/'

gatk="/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk"

# set number of jobs based on number of cores available and number of threads per job
export JOBS_PER_NODE=$(( $SLURM_CPUS_ON_NODE / $SLURM_CPUS_PER_TASK ))

#echo $SLURM_JOB_NODELIST |sed s/\,/\\n/g > hostfile

# Assuming this is being run in the same directory as ChromNames.txt and the shell script. 

for samp in $(cat ../bam-sets/bqsr-bams-set-00)
do
	id=$(echo $samp | sed 's/_BQSR.bam//g')
	parallel \
	--jobs $JOBS_PER_NODE \
	--joblog ./extractVars-out/extract-vars-${id}.log \
	--resume --progress \
	-a ./ChromNames.txt ./GATK-Extract-VarsByChrom.sh {} $id $SLURM_CPUS_PER_TASK
done
