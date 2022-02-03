#!/bin/bash
#SBATCH --output=./err-out/format.msmc_%A_%a.out
#SBATCH --error=./err-out/format.msmc_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --array=1-11
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=fc_fishes
#SBATCH -p savio3

# This script will loop through the per-scaffold vcfs for high-coverage samples, 
# converting vcfs to the input format ('multihetsep') for MSMC. 
# This is done in parallel such that each sample is being processed simultaneously

module load bcftools
module load samtools
module load python/3.7

id=$SLURM_ARRAY_TASK_ID
samp=$(sed "${id}q;d" msmc-samps.txt)

echo "Preparing MSMC input for $samp"

MSMCTOOLS=/global/scratch/users/austinhpatton/software/msmc-tools
maskDir=/global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/masks
msmcIn=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/MSMC/Input
msmcSnps=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/SamtoolsGenotypes

# Loop through each linkage group
for scaff in $(cat ChromNames.txt)
do 
	echo "Converting $scaff to multihetsep"
	
	# Set the path/filename to which youll be saving the input to a variable
	multihet=$msmcIn/$samp-$scaff-msmc-multihetsep.txt

	# Store the path to the VCF for this sample and associated mask-file to variables
	maskInd=$maskDir/IndMask-$samp-$scaff.bed.gz
	maskGenome="$maskDir/mask_50_50-$scaff.mask.bed.gz"
	vcf="$msmcSnps/$samp-$scaff.vcf.gz"

	# Now, convert the vcf to multihetsep msmc input using both an individual mask and genome-wide mappability mask
	$MSMCTOOLS/generate_multihetsep.py --mask=$maskInd --mask=$maskGenome $vcf > $multihet

	echo "Finished with $scaff. Moving on." 
done
echo "Data formatted for $samp"
