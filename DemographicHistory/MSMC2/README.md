## MSMC2 Presp & Analysis Scripts

This directory holds a number of scripts to be run in order. A description of the steps are outlined below. 

### 1) Generate a mappability mask where you chop up the reference into "reads", map those back to the reference, and only retain sites (non-intuitively, these are the sites that are 'masked') that can be reliably called using short reads. This involves running "Make-SNPable-Mask.slurm" and then "Make-Mappability-Mask.slurm" using the mask file the first script produces. 
### 2) Then, call sites per individual (and per scaffold) using samtools (msmc-callSites-1.slurm), and 
### 4) convert to the msmc input file (msmc-formatInput-2.slurm). 
#### - These scripts are written with an emphasis on parallelization, so these two scripts run as an array job on slurm schedulers, where each job within the array is processing a different sample, and then within the scripts you loop through scaffolds/chromosomes. 
### 5) After you've got the input files made, you run "msmc-runPerInd-3.slurm" which does the full analysis as an array job, again with individuals being the unit of parallelization, and providing all input files you made. 
### 6) Lastly, you run "msmc-bootPerInd-4.slurm", again parallelized by individual, which does the whole bootstrapping process - resampling scaffolds, etc. 
