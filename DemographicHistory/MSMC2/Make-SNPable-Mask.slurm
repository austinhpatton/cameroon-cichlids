#!/bin/bash
#SBATCH --output=snpable_%a.out
#SBATCH --error=snpable_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=50:00:00
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH -p savio3
#SBATCH --account=co_fishes

module load python
module load samtools
module load bcftools
module load bwa 
module load java

REF='/global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fa'
MASK='/global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.mask.fa'
IDX='/global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fai'
OUT='/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/SamtoolsGenotypes/'
bamDir='/global/scratch/users/austinhpatton/cichlids/cameroon/preProcessing/mappedBams/'

gatk="/global/scratch/users/austinhpatton/software/gatk-4.1.8.1/gatk"

cd /global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/

seqbility="/global/scratch/users/austinhpatton/software/seqbility-20091110"

# Extract overlapping kmer subsequences as reads
$seqbility/splitfa $REF 50 | split -l 20000000

# This splits the "reads into a bunch of files. Combine into a single "read" file and clean up
cat x* > split-fa-reads.fa
rm x* 

# Now map these "Reads" back to the genome
bwa aln -t 30 -R 1000000 -O 3 -E 3 $REF split-fa-reads.fa > split-fa-reads.sai  

# Convert to sam
bwa samse -f split-fa-reads.sam $REF split-fa-reads.sai split-fa-reads.fa

# Generate the raw mask from the unsorted sam index 
gzip split-fa-reads.sam
gzip -dc split-fa-reads.sam.gz | $seqbility/gen_raw_mask.pl > rawMask_50.fa  

# and finally generate the full mask
$seqbility/gen_mask -l 50 -r 0.5 rawMask_50.fa > mask_50_50.fa  

# change all bases that are not mapping uniquely and with some degree of mismatch (c!=3) to lowercase
$seqbility/apply_mask_s mask_50_50.fa $REF > $MASK  
