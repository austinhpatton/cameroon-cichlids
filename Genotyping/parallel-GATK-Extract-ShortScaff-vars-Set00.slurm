#!/bin/bash
#SBATCH --output=ExtractScaffVars_set-00_%a.out
#SBATCH --error=ExtractScaffVars_set-00_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --exclusive
#SBATCH --account=ac_fishes
#SBATCH -p savio3

# Interpretation/functionality is the same as in parallel-GATK-Extract-LongScaffs-vars-Set00.slurm.
# Difference is, here we do not loop through samples, running scaffolds in parallel. Here, we analyze 
# all short scaffolds together, running samples in parallel. 
# Again, this is really to circumvent the 72 hour limit on savio compute times. 

# Load the default version of GNU parallel.
module load gnu-parallel
module load java 

REF='/global/scratch/austinhpatton/cichlids/Oreochromis-Reference/Oreochromis.fna'
INTERIM='/global/scratch/austinhpatton/cichlids/cameroon/Interim-BQSR/'
bamDir='/global/scratch/austinhpatton/cichlids/cameroon/preProcessing/mappedBams/'
gatkOut='/global/scratch/austinhpatton/cichlids/cameroon/gatk-output/bqsr-prep/'

gatk="/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk"

# set number of jobs based on number of cores available and number of threads per job
export JOBS_PER_NODE=$(( $SLURM_CPUS_ON_NODE / $SLURM_CPUS_PER_TASK ))

parallel \
	--jobs $JOBS_PER_NODE \
	--joblog extractVars-out/extract-vars-markdup-00-short.log \
	--resume --progress \
	-a markdup-set-00 ./GATK-Extract-VarsByScaff.sh {} $SLURM_CPUS_PER_TASK
