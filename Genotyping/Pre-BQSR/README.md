This folder contains all scripts used for the first round of genotyping as well as base quality score recalibration. 

Scripts are ran as follows:
1) parallel-GATK-Extract-LongScaffs-vars-Set00.slurm is used to submit GATK-Extract-VarsByLongScaff.sh using gnu-parallel
2) parallel-GATK-Extract-ShortScaffs-vars-Set00.slurm is used to submit GATK-Extract-VarsByShortScaff.sh using gnu-parallel
3) GatherGVCFs.slurm is used to submit GatherGVCFs.sh, in all generating per-sample gvcfs, 20 samples in parallel
4) GATK-GenotypeGVCFs.slurm is not parallelized (boo joint genotyper), but conducts joint genotyping, pull out indels and snps, and hard filters them
5) GATK-RecalModel-Parallel.slurm runs GATK-BaseRecalibrator.sh to build BQSR models for all samples in parallel
6) GATK-ApplyBQSR-Parallel.slurm applies that BQSR model using GATK-ApplyBQSR.sh to spit out a recalibrated bam for each sample
