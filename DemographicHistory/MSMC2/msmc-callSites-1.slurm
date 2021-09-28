#!/bin/bash
#SBATCH --output=map.msmc_%A_%a.out
#SBATCH --error=map.msmc_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --array=1-11
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=co_rosalind
#SBATCH -p savio3

# This script will loop through the high-coverage samples, mapping reads
# using bcftools in preparation for use with MSMC. 

module load bcftools
module load samtools
module load python/3.7

id=$SLURM_ARRAY_TASK_ID
samp=$(sed "${id}q;d" msmc-samps.txt)

echo "Processing data for ${samp}"

MSMCTOOLS=/global/scratch/users/austinhpatton/software/msmc-tools
bam=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/RecalBams/${samp}_BQSR.bam
maskOut=/global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/masks
msmcOut=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/SamtoolsGenotypes
REF=/global/scratch/users/austinhpatton/cichlids/Oreochromis-Reference-2018/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fa

if [ -f "${bam}.bai" ]
        then
                echo "Bamfile index exists"
        else
                echo "Creating index"
                samtools index $bam ${bam}.bai
fi

#Now, we need to call snps for each linkage group
for scaff in $(cat ChromNames.txt)
do 
	echo "Beginning scaffold $scaff"

	# Calculate mean coverage for use in msmc-tools' bamCaller.py script
	cov=$(samtools depth -r $scaff $bam | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.')
	echo "Mean coverage for $scaff is $cov"

	# Now, generate the VCF for this sample and an associated mask-file
	mask=$maskOut/IndMask-${samp}-${scaff}.bed.gz
	vcf="${msmcOut}/${samp}-${scaff}.vcf"

	# Now, call SNPs using samtools, excluding low coverage (<10X) sites. 
	bcftools mpileup -Ou -r $scaff --threads 30 -f $REF $bam | bcftools call -c --threads 30 -V indels | $MSMCTOOLS/bamCaller.py $cov $mask | gzip -c > ${vcf}.gz	
	#samtools mpileup -q 20 -Q 20 -C 50 -u -r $scaff -f $REF $bam | bcftools call -c -V indels | bcftools view -i 'INFO/DP>9' | $MSMCTOOLS/bamCaller.py $cov $mask | gzip -c > ${vcf}.gz

	echo "Finished with ${scaff}. Moving on." 
done
echo "Genotype calling completed for $samp"
