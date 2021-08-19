#!/bin/bash
#SBATCH --output=GatherVCFs_%a.out
#SBATCH --error=GatherVCFs_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=5000:00:00
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --account=co_fishes
#SBATCH -p savio3
###SBATCH --mem=96G

# Load the default version of GNU parallel.
module load gnu-parallel
module load java 

# This script run GATK GatherVCFs, 20 samples at a time using gnu-parallel. 
parallel --jobs 20 -a BMsamps.txt ./GatherGVCFs.sh {}

