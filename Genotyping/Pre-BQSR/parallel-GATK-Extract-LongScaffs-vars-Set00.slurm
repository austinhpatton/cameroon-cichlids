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

# This is one of several parallel submission scripts used to submit "GATK-Extract-LongScaffs.sh". 
# Several of these scripts were used, each on "markdup-set-XX", each including a different set of samples. 
# In an unconstrained cluster, this could be used on the full set of samples, but I needed to break these 
# up so that they would finish within 72 hours and also to navigate computational (i.e. memory) limitations. 
# In short, this script loops through samples, calculating genotype likelihoods for the long-scaffolds 
# in "ChromNames.txt" in parallel. Embarrassingly parallel.

# Load the default version of GNU parallel.
module load gnu-parallel
module load java

# set number of jobs based on number of cores available and number of threads per job
export JOBS_PER_NODE=$(( $SLURM_CPUS_ON_NODE / $SLURM_CPUS_PER_TASK ))

#echo $SLURM_JOB_NODELIST |sed s/\,/\\n/g > hostfile

# markdup-set-00 contains a reduced list of "markdup.bam"s to be run. 
# "ChromNames.txt" is a list of the 22 long, roughly chromosome length scaffolds in the reference. 
# By feeding these to parallels, we conduct a scatter-call operation with GATK.
# In other words, we loop through samples, calculating genotype likelihoods for long scaffolds in parallel. 
# Together, this expedites the genotyping process. 
for samp in $(cat markdup-set-00)
do
  	id=$(echo $samp | sed 's/_markdup.bam//g')
        parallel \
        --jobs $JOBS_PER_NODE \
        --joblog extractVars-out/extract-vars-${id}.log \
        --resume --progress \
        -a ChromNames.txt ./GATK-Extract-VarsByChrom.sh {} $id $SLURM_CPUS_PER_TASK
done



