# Divergent selection using CNNs
Scripts pertaining to the simulation of loci under divergent selection, sweeps, or neutrality using the demographic histories estimated for focal species. 

### Note that many of these scripts are modified from those found here: 
https://github.com/alexnater/midas-genomics/tree/master/divergent_selection 

### Scripts include the following:
1) generate_tbs.py - this is used by generateSimLaunchScrips_msms.py to draw random samples of parameters in the simulations so as to introduce stochasticity into the coalescent simulations
2) generate_sims_msms.sh this calls generateSimLaunchScripts_msms.py given a set of parameters (i.e. those above) to generate a set of scripts specifying simulations with msms
3) generateSimLaunchScripts_msms.py - this is called by the above scripts, and generates a set of coalescent simulations to be used in training. 
