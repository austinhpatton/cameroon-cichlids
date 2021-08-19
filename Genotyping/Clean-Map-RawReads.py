#!/usr/bin/
# Full pipeline for everything from read merging, 
# trimming, mapping to a reference, sorting resultant 
# bam files, and marking/removing PCR duplicates.


#bwa reference-based assembly illumina Paired End WGS data

from os.path import join as jp
from os.path import abspath
import os
import sys
import argparse

#sys.argv
parser = argparse.ArgumentParser()
parser.add_argument('-s', "--samples", help="Samples.txt file with sample ID.", required=True)
parser.add_argument('-r', "--rawdata", help="Path to raw fastq data.", required=True)
parser.add_argument('-b', "--bwaindex", help="Path to bwa index file.", required=True)
args = parser.parse_args()

VERBOSE=True

#Function definitions:
def log(txt, out):
    if VERBOSE:
        print(txt)
    out.write(txt+'\n')
    out.flush()

## Read in samples and put them in a list:
samples = []
for l in open(args.samples):
    if len(l) > 1:
        samples.append(l.split('_R')[0].replace('_R1_001.fastq.gz', '').strip())

# Setup folders and paths variables:
resultsDir = "/global/scratch/austinhpatton/cichlids/cameroon/preProcessing/cleanedReads"
bamFolder = "/global/scratch/austinhpatton/cichlids/cameroon/preProcessing/mappedBams"
PBS_scripts = "/global/scratch/austinhpatton/cichlids/cameroon/preProcessing/CleaningScripts"
rawdataDir = abspath(args.rawdata)
bwaIndex = abspath(args.bwaindex)
picard="/global/home/users/austinhpatton/software/java-1.8/bin/java -jar /global/home/users/austinhpatton/software/picard.jar "
#bwa="/home/matthew.lawrance/bwa-0.7.17/bwa "
#flash2="/global/home/users/austinhpatton/software/flash2/flash2 "
#sickle="/global/home/users/austinhpatton/software/sickle/sickle "

os.system('mkdir -p %s' % resultsDir)
os.system('mkdir -p %s' % bamFolder)
os.system('mkdir -p %s' % PBS_scripts)

##### Run pipeline ###
for sample in samples:
    # Set up files:
    logFile = jp(resultsDir, sample + '_cleaning.log')
    logCommands = open(jp(PBS_scripts, sample + '_cleaning_commands.sh'), 'w')
# Merge reads using flash2
    cmd = ' '.join(['/global/home/users/austinhpatton/software/flash2/flash2 --max-overlap 300 --allow-outies --threads 7', ' -d ', resultsDir, ' -o ', jp(sample + '_flash'),
                    jp(rawdataDir, sample + '_R1_001.fastq.gz'), jp(rawdataDir, sample + '_R2_001.fastq.gz'),
                    '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    os.system(cmd)
# Trim reads using sickle
    cmd = ' '.join(['/global/home/users/austinhpatton/software/sickle/sickle pe --length-threshold 20 --qual-threshold 25 --qual-type sanger -f', jp(resultsDir, sample + '_flash.notCombined_1.fastq'),
                    '-r', jp(resultsDir, sample + '_flash.notCombined_2.fastq'),
                    '--output-pe1', jp(resultsDir, sample + '_sickle_PE1.fastq'),
                    '--output-pe2', jp(resultsDir, sample + '_sickle_PE2.fastq'),
                    '--output-single', jp(resultsDir, sample + '_sickle_SE.fastq'), '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    os.system(cmd) 
# Combine non-merged, single-end read files:
    cmd = ' '.join(['cat', jp(resultsDir, sample + '_sickle_SE.fastq'), jp(resultsDir, sample + '_flash.extendedFrags.fastq'),
                    '>', jp(resultsDir, sample + "_cleaned_SE.fastq")])
    log(cmd, logCommands)
    os.system(cmd)
# Rename PE and SE files:
    cmd = ' '.join(['mv', jp(resultsDir, sample + "_sickle_PE1.fastq"), jp(resultsDir, sample + "_cleaned_PE1.fastq")])
    log(cmd, logCommands)
    os.system(cmd)
    cmd = ' '.join(['mv', jp(resultsDir, sample + "_sickle_PE2.fastq"), jp(resultsDir, sample + "_cleaned_PE2.fastq")])
    log(cmd, logCommands)
    os.system(cmd)
# Clean up intermediary files:
    cmd = ' '.join(['rm', jp(resultsDir, sample + "_sickle*"), jp(resultsDir, sample + "_flash.extendedFrags.fastq")])
    log(cmd, logCommands)
    os.system(cmd)
# Compress cleaned files:
    cmd = ' '.join(['gzip', jp(resultsDir, sample + '*.fastq')])
    log(cmd, logCommands)
    os.system(cmd)
# Run BWA to map PE samples to reference genome
# -t number of threads -R read group header
    logFile = jp(resultsDir, sample + '_mapping.log')
    cmd = ' '.join(["/global/home/users/austinhpatton/software/bwa/bwa mem -t 50 -R '@RG\\tID:bwa\\tSM:" + sample + "\\tPL:ILLUMINA'",
                    bwaIndex, jp(resultsDir, sample + "_cleaned_PE1.fastq.gz"),
                    jp(resultsDir, sample + "_cleaned_PE2.fastq.gz"), "| /global/home/users/austinhpatton/software/samtools-1.11/samtools view -bS -@ 50 -o", jp(bamFolder, sample + "PE.bam"),
                    "2>", logFile])
    log(cmd, logCommands)
    os.system(cmd)
# sort PE bam file; -@ number of threads
    cmd = ' '.join(["/global/home/users/austinhpatton/software/samtools-1.11/samtools sort -o", jp(bamFolder, sample) + "sortedPE.bam", ' -@ 60', jp(bamFolder, sample + "PE.bam")])
    log(cmd, logCommands)
    os.system(cmd)
# Run BWA to map SE samples to reference genome
    cmd = ' '.join(["/global/home/users/austinhpatton/software/bwa/bwa mem -t 50 -R '@RG\\tID:bwa\\tSM:" + sample + "\\tPL:ILLUMINA'",
                    bwaIndex, jp(resultsDir, sample + "_cleaned_SE.fastq.gz"), "| /global/home/users/austinhpatton/software/samtools-1.11/samtools view -bS -@ 50 -o", jp(bamFolder, sample + "SE.bam"),
                    "2>>", logFile])
    log(cmd, logCommands)
    os.system(cmd)
# sort SE bam file; -@ number of threads
    cmd = ' '.join(["/global/home/users/austinhpatton/software/samtools-1.11/samtools sort -o", jp(bamFolder, sample) + "sortedSE.bam", ' -@ 60', jp(bamFolder, sample + "SE.bam")])
    log(cmd, logCommands)
    os.system(cmd)    
# merge
    cmd = ' '.join(["/global/home/users/austinhpatton/software/samtools-1.11/samtools merge -c", jp(bamFolder, sample + ".bam"), jp(bamFolder, sample + "sortedPE.bam"), jp(bamFolder, sample + "sortedSE.bam")])
    log(cmd, logCommands)
    os.system(cmd)
# sort bam file; -@ number of threads
    cmd = ' '.join(["/global/home/users/austinhpatton/software/samtools-1.11/samtools sort -o", jp(bamFolder, sample) + "sorted.bam", ' -@ 60', jp(bamFolder, sample + ".bam")])
    log(cmd, logCommands)
    os.system(cmd)
# Mark and remove PCR duplicates
    cmd = ' '.join([picard + "MarkDuplicates INPUT=" + jp(bamFolder, sample + "sorted.bam"), ' OUTPUT=' + jp(bamFolder, sample + "_markdup.bam"),
                    ' METRICS_FILE=' + jp(bamFolder, sample + ".metrics"), ' REMOVE_DUPLICATES=true ',
                    ' ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT', '>>', logFile, '2>&1'])
    log(cmd, logCommands)
    os.system(cmd)
# Index bam file:
    cmd = ' '.join(["/global/home/users/austinhpatton/software/samtools-1.11/samtools index", jp(bamFolder, sample) + "_markdup.bam"])
    log(cmd, logCommands)
    os.system(cmd)   
    logCommands.close()

