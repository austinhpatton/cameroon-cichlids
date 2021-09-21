#!/bin/bash
#SBATCH --output=./err-out/FitRecal_%a.out
#SBATCH --error=./err-out/FitRecal_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --account=fc_fishes
#SBATCH -p savio3

# Load the default version of GNU parallel.
module load gnu-parallel
module load java 

# This script fits the base score quality recalibration tool to each sample, in parallel. 
parallel --jobs 10 -a BM-Ckott-samps.txt ./GATK-BaseRecalibrator.sh {}

