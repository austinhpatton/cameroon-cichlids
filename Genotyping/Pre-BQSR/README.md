This folder contains all scripts used for the first round of genotyping as well as base quality score recalibration. 

Note that many of these scripts are parallelized in a manner to circumvent the wall time limit (72 hrs) imposed by the savio cluster. 
Without such a limit, rather than running, for instance, parallel-GATK-Extract-ShortScaffs-vars-Set00.slurm on a set of 4-5 samples, this could
be done using a complete set of samples. 

There are unfortunately more steps doing it this way, although it is substantially more efficient!

Scripts are ran as follows:
1) parallel-GATK-Extract-LongScaffs-vars-Set00.slurm is used to submit GATK-Extract-VarsByLongScaff.sh using gnu-parallel
2) parallel-GATK-Extract-ShortScaffs-vars-Set00.slurm is used to submit GATK-Extract-VarsByShortScaff.sh using gnu-parallel
3) DBImport-ByScaff.slurm generates a genomic database per scaffold/chromosome for more efficient joint genotyping 
4) GATK-GenotypePerScaffGVCFs.slurm conducts joint genotyping per scaffold/chromosome
5) GatherJointVCFs.slurm combined these scatter-called joint vcfs into one
6) GATK-ExtractFilter-Variants.slurm will take the combined, joint genotyped VCF and pull out indels and snps, and hard filter them
7) GATK-RecalModel-Parallel.slurm runs GATK-BaseRecalibrator.sh to build BQSR models for all samples in parallel
8) GATK-ApplyBQSR-Parallel.slurm applies that BQSR model using GATK-ApplyBQSR.sh to spit out a recalibrated bam for each sample

Then, you repeat 1-6 and voila, you've got a set of variants based on recalibrated base quality scores!
