#!/bin/bash
#SBATCH --output=FitRecal_%a.out
#SBATCH --error=FitRecal_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --account=co_fishes
#SBATCH -p savio3
###SBATCH --mem=96G

# Load the default version of GNU parallel.
module load gnu-parallel
module load java 

# This script fits the base score quality recalibration tool to each sample, in parallel. 
parallel --jobs 10 -a BMsamps.txt ./GATK-BaseRecalibrator.sh {}

