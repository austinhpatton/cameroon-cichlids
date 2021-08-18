#!/bin/bash
#SBATCH --output=preProcess_%J.out
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --array=1-26
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio3
#SBATCH --account=ac_fishes

module load python
module load samtools
module load bcftools
module load bwa 

# The raw sequences for each species are in separate directories, so we'll take advantage of that and just do the preprocessing by species. 
# In all that will be 26 independent sets of runs submitted as an array, and the resulting files will all go to a single directory.
id=$SLURM_ARRAY_TASK_ID

# Extract the target directory holding raw reads using that array index id:
dir=$(sed -n ${id}p preprocess-interim/CameroonCichs-RawSeqDirs.txt)

# Now, pull out species name into a variable to use in file output
# Each directory follows the same pattern - 
# /global/scratch/austinhpatton/cichlids/cameroon/raw-seqs/SPECIES_NAME
# So, extract the species ID using cut - using '/' as the delimiter, species is the 8th field
spp=$(echo $dir | cut -d'/' -f8)

ls ${dir}/*fastq.gz | cut -d'/' -f9 > RawReadLists/${spp}-raw-reads.txt
echo $spp

# Note that below, "${spp}-raw-reads.txt" is a file containing a list of file names for raw read for each species. 
python Clean-Map-RawReads.py -s RawReadLists/${spp}-raw-reads.txt -r ${dir}/ -b /global/scratch/austinhpatton/cichlids/Oreochromis-Reference/Oreochromis.fna 
