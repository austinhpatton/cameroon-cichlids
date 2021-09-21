This folder contains all scripts used for the first round of genotyping as well as base quality score recalibration. 

Note that many of these scripts are parallelized in a manner to circumvent the wall time limit (72 hrs) imposed by the savio cluster. 
Without such a limit, rather than running, for instance, parallel-GATK-Extract-ShortScaffs-vars-Set00.slurm on a set of 4-5 samples, this could
be done using a complete set of samples. 

Scripts are ran as follows:
1) parallel-GATK-Extract-LongScaffs-vars-Set00.slurm is used to submit GATK-Extract-VarsByLongScaff.sh using gnu-parallel
2) parallel-GATK-Extract-ShortScaffs-vars-Set00.slurm is used to submit GATK-Extract-VarsByShortScaff.sh using gnu-parallel
3) DBImport-ByScaff.slurm generates a genomic database per scaffold/chromosome for more efficient joint genotyping 
4) GATK-GenotypeGVCFs.slurm is not parallelized (boo joint genotyper), but conducts joint genotyping, pull out indels and snps, and hard filters them
5) GATK-RecalModel-Parallel.slurm runs GATK-BaseRecalibrator.sh to build BQSR models for all samples in parallel
6) GATK-ApplyBQSR-Parallel.slurm applies that BQSR model using GATK-ApplyBQSR.sh to spit out a recalibrated bam for each sample
