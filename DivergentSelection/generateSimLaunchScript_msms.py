#!/usr/bin/env python3

import sys, os, math, statistics

"""
Prints a bash script to stdout for running training and testing simulations with msms.
"""


def read_params(infile, pops):
    with open(infile, 'r') as inhandle:
        # read header line:
        names = inhandle.readline().split()[1:]
        values = map(float, inhandle.readline().split()[1:])
    
    # replace populations labels by numbers:
    if len(pops) != 3:
        raise Exception("ERROR: Provide list with 3 population labels!")
    for popi in range(len(pops)):
        names = [s.replace(pops[popi], "pop" + str(popi + 1), 1) for s in names]
    
    return dict(zip(names, values))

# a function to get the popsize history for each species
def get_popHist(popSizeFileName):
    sizeChanges = []
    with open(popSizeFileName) as popSizeFile:
        first = True
        for line in popSizeFile:
            year, ne = line.strip().split()
            year, ne = float(year), float(ne)
            if first:
                first = False
                ne0 = ne
                prevSizeRatio = 1.0
            else:
                t, sizeRatio = year*gensPerYear/(4*ne0), ne/ne0
                if abs(sizeRatio - prevSizeRatio) > 1e-9:
                    sizeChanges.append((t, sizeRatio, ne))
                    prevSizeRatio = sizeRatio
    return sizeChanges, ne0
    
# python script to generate random parameter values from specified distributions.
tbs_gen = "/data/scc3/anater/output/diffstats/msms_mom/generate_tbs.py"

trainingOutDir = sys.argv[1]	# output directory for training data
testOutDir = sys.argv[2]	# output directory for testing data
outPrefix = sys.argv[3]	# prefix for output files
trainingSampleNumber = int(sys.argv[4])	# number of simulation replicates for training data
testSampleNumber = int(sys.argv[5])	# number of simulation replicates for testing data
sampleSizes = [ int(x) for x in sys.argv[6].split(',') ]	# number of individuals in each population sample
numSites = int(sys.argv[7])	# total number of sites in our simulated window
rhoOverTheta = float(sys.argv[8])	# ratio of recombination to mutation rate
### THIS IS WHAT I CHANGE TO INCLUDE THE MSMC (OR WHATEVER) POPULATION SIZE HISTORY
# RATHER THAN SIMULATING WITH MIGRATION, ETC, JUST USE THE POPSIZEHISTORY
popsize_1 = sys.argv[9] # pop size history file for spp 1
popsize_2 = sys.argv[10] # pop size history file for spp 2
# param_file = sys.argv[9]	# file with maximum likelihood parameter estimates from fastsimcoal # CHANGE TO JUST NE, GENS BEFORE PRESENT
###
pops = sys.argv[11].split(',')	# list of population abbreviations for later use # MAYBE NOT?
tdiv = ses.argv[12] # time at which species diverged

# hard-coded settings:
u = 3.5e-9	# per-site mutation rate
freqLow = 1e-5	# minimum initial frequency of selected allele
freqHigh = 0.01	# maximum initial frequency of selected allele
sLow =0.01	# minium selection coefficient 
sHigh = 0.1	# maximum selection coefficient 
posLow = 0.4	# minimum of position of selected site
posHigh = 0.6	# maximum of position of selected site
minsweepfreq = 0.8	# minimum frequency of selected allele. AFTER SELECTION?
gensPerYear = 0.75  # generations per year (i.e. generation time is 1.5 years)

sampleSize = sum(sampleSizes)
thetaMean = 4 * N0 * u * numSites
rhoMean = thetaMean * rhoOverTheta
rhoMax = 3 * rhoMean
thetaLow = (2 * thetaMean) / 11.0
thetaHigh = 10 * thetaLow
alphaLow = 2 * N0 * sLow
alphaHigh = 2 * N0 * sHigh

# read in population size histories
sizeChanges_1, Ne1_curr = get_pophist(popsize_1)
sizeChanges_2, Ne2_curr = get_pophist(popsize_2)
# calculate ancestral population size as the mean of the oldest pop sizes for each species
Ne_anc = mean(sizeChanges_1[16][2], sizeChanges_2[16][2])
# calculate population-wide Ne at present as the sum of the oldest pop sizes for each species
N0 = sizeChanges_1[0][2] + sizeChanges_2[0][2]



params = read_params(param_file, pops)
tsplit_cl = params["TDIV_cl"] / (4*N0)
# tadmix = params["TADMIX"] / (4*N0)
tcol = params["TDIV_col"] / (4*N0)
tsplit_gl = params["TDIV_gl"] / (4*N0)
tchange = params["TCHANGE_pop1"] / (4*N0)

sharedSelStr = "-Sc 0 1 tbs tbs 0 -Sc 0 2 tbs tbs 0 -SI {} 2 0 0 tbs tbs -Sp tbs -SFC -oTrace -oOC".format(tadmix)
indepSelStr = "-Sc 0 1 tbs tbs 0 -Sc 0 2 tbs tbs 0 -SI {} 2 0 0 tbs tbs -Sp tbs -SFC -oTrace -oOC".format(tsplit_cl)
div1SelStr = "-Sc 0 1 tbs tbs 0 -Sc 0 2 tbs tbs 0 -SI {} 2 0 0 tbs tbs -Sp tbs -SFC -oTrace -oOC".format(tsplit_cl)
div2SelStr = "-Sc 0 1 tbs tbs 0 -Sc 0 2 tbs tbs 0 -SI {} 2 0 0 tbs tbs -Sp tbs -SFC -oTrace -oOC".format(tsplit_cl)
demogStr1 = "-I 2 {} {} -n 1 {} -n 2 {}".format(*sampleSizes, Ne1_curr, Ne2_curr)
demogStr2 = "-ej {} 4 3 -es {} 3 {} -ej {} 5 1 -es {} 3 {} -ej {} 6 2".format(tsplit_cl, tadmix, 1-params["MIGRANTS_pop1"], tadmix, tadmix, 1-params["MIGRANTS_pop2"], tadmix)
demogStr3 = "-es {} 3 {} -ej {} 7 2 -ej {} 3 1".format(tcol, 1-params["CONT_pop2"], tcol, tcol)
demogStr4 = "-ej {} 2 1 -en {} 1 {}".format(tsplit_gl, tchange, Ne_anc)

neutFilterStr = 'awk \'/^ms/{next;} /^0x/{next;} 1{print}\''
sharedFilterStr = 'awk -v minfreq=%f \'BEGIN{validsim=0} /^ms/{next} /^\/\//{output=$0; validsim=1; readseqs=0; freq1=0; freq2=0; next} !validsim{next} NF==9{freq1=$7; freq2=$9; next} /^segsites/{output=output"\\n"$0; if(freq1<minfreq || freq2<minfreq) validsim=0; next} /^positions/{output=output"\\n"$0; readseqs=1; next} readseqs && /^[01]/{output=output"\\n"$0; next} /^OriginCount/{if(validsim && readseqs) print output"\\n"; next}\'' %(minsweepfreq)
indepFilterStr = 'awk -v minfreq=%f \'BEGIN{validsim=0} /^ms/{next} /^\/\//{output=$0; validsim=1; readseqs=0; freq1=0; freq2=0; next} !validsim{next} NF==9{freq1=$7; freq2=$9; next} /^segsites/{output=output"\\n"$0; if(freq1<minfreq || freq2<minfreq) validsim=0; next} /^positions/{output=output"\\n"$0; readseqs=1; next} readseqs && /^[01]/{output=output"\\n"$0; next} /^OriginCount/{if(validsim && readseqs) print output"\\n"; next}\'' %(minsweepfreq)
div1FilterStr = 'awk -v minfreq=%f \'BEGIN{validsim=0} /^ms/{next} /^\/\//{output=$0; validsim=1; readseqs=0; freq1=0; freq2=0; next} !validsim{next} NF==9{freq1=$7; freq2=$9; next} /^segsites/{output=output"\\n"$0; if(freq1<minfreq) validsim=0; next} /^positions/{output=output"\\n"$0; readseqs=1; next} readseqs && /^[01]/{output=output"\\n"$0; next} /^OriginCount/{if(validsim && readseqs) print output"\\n"; next}\'' %(minsweepfreq)
div2FilterStr = 'awk -v minfreq=%f \'BEGIN{validsim=0} /^ms/{next} /^\/\//{output=$0; validsim=1; readseqs=0; freq1=0; freq2=0; next} !validsim{next} NF==9{freq1=$7; freq2=$9; next} /^segsites/{output=output"\\n"$0; if(freq2<minfreq) validsim=0; next} /^positions/{output=output"\\n"$0; readseqs=1; next} readseqs && /^[01]/{output=output"\\n"$0; next} /^OriginCount/{if(validsim && readseqs) print output"\\n"; next}\'' %(minsweepfreq)

shared_mult = 1.05
indep_mult = 1.05
div_mult = 1.05

print("#!/bin/bash")
for sampleNumber, outDir, simTitle in [(trainingSampleNumber, trainingOutDir, "training data"), (testSampleNumber, testOutDir, "test data")]:
    print("\n# generating %s\n" %(simTitle))
    neutMSMSCmd = "python %s %d %d - -t %f %f -r %f %f | msms %d %d -N %d -t tbs -r tbs %d" %(tbs_gen, 2, sampleNumber, thetaLow, thetaHigh, rhoMean, rhoMax, sampleSize, sampleNumber, N0, numSites)
    print("{} {} {} {} {} | {} > {}/{}.neut.msOut".format(neutMSMSCmd, demogStr1, demogStr2, demogStr3, demogStr4, neutFilterStr, outDir, outPrefix))  
    
    sharedMSMSCmd = "python %s %d %d - -t %f %f -r %f %f -a %f %f %f %f -f %f %f %f %f -p %f %f --log | msms %d %d -N %d -t tbs -r tbs %d" %(tbs_gen, 2, sampleNumber * shared_mult, thetaLow, thetaHigh, rhoMean, rhoMax, alphaLow, alphaHigh, alphaLow, alphaHigh, freqLow, freqHigh, 0., 0., posLow, posHigh, sampleSize, sampleNumber * shared_mult, N0, numSites)
    indepMSMSCmd = "python %s %d %d - -t %f %f -r %f %f -a %f %f %f %f -f %f %f %f %f -p %f %f --log | msms %d %d -N %d -t tbs -r tbs %d" %(tbs_gen, 2, sampleNumber * indep_mult, thetaLow, thetaHigh, rhoMean, rhoMax, alphaLow, alphaHigh, alphaLow, alphaHigh, freqLow, freqHigh, freqLow, freqHigh, posLow, posHigh, sampleSize, sampleNumber * indep_mult, N0, numSites)
    div1MSMSCmd = "python %s %d %d - -t %f %f -r %f %f -a %f %f %f %f -f %f %f -p %f %f --log | msms %d %d -N %d -t tbs -r tbs %d" %(tbs_gen, 2, sampleNumber * div_mult, thetaLow, thetaHigh, rhoMean, rhoMax, alphaLow, alphaHigh, -alphaLow, -alphaHigh, freqLow, freqHigh, posLow, posHigh, sampleSize, sampleNumber * div_mult, N0, numSites)
    div2MSMSCmd = "python %s %d %d - -t %f %f -r %f %f -a %f %f %f %f -f %f %f -p %f %f --log | msms %d %d -N %d -t tbs -r tbs %d" %(tbs_gen, 2, sampleNumber * div_mult, thetaLow, thetaHigh, rhoMean, rhoMax, -alphaLow, -alphaHigh, alphaLow, alphaHigh, freqLow, freqHigh, posLow, posHigh, sampleSize, sampleNumber * div_mult, N0, numSites)
    print("")
    print("{} {} {} {} {} {} | {} > {}/{}.shared.msOut".format(sharedMSMSCmd, demogStr1, demogStr2, demogStr3, demogStr4, sharedSelStr, sharedFilterStr, outDir, outPrefix))
    print("{} {} {} {} {} {} | {} > {}/{}.divergent1.msOut".format(div1MSMSCmd, demogStr1, demogStr2, demogStr3, demogStr4, div1SelStr, div1FilterStr, outDir, outPrefix))
    print("{} {} {} {} {} {} | {} > {}/{}.divergent2.msOut".format(div2MSMSCmd, demogStr1, demogStr2, demogStr3, demogStr4, div2SelStr, div2FilterStr, outDir, outPrefix))
    print("")

