#!/bin/bash
#SBATCH --output=bqsr_%a.out
#SBATCH --error=bsqr_%a.err
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

# This script submits the shell script that applies the base quality score recalibration to each samples bam file (listed in BMsamps.txt), in parallel. 

parallel --jobs 10 -a BMsamps.txt ./GATK-ApplyBQSR.sh {}

