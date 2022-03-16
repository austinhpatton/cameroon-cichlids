# This script takes in the results from Dinvestigate and spits out a nice 
# little plot of various D-statistics from the sliding windows
setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/cichlids/cameroon/Onil_UMD/Dsuite/Results/100Snp_Slide25/FocalTrios/')
library(ggplot2)
library(inlmisc)
library(scales)
library(ggallin)
library(cowplot)

###########################################################################################
# For each species, we want to sample the window who's midpoint has the highest  
# f_dm statistics per 5 kb

# So, we want to focus on the strongest cases for introgression for each of our focal species.
# Specifically, we want to look at genome-wide signatures of introgression either 
# 1) in sympatry (e.g. between species within Barombi Mbo), or
# 2) involving allopatric species (e.g. between Barombi Mbo and S. galileus 
# from either Meme/Mungo or Cross River)

# Using df as our statistic of choice, the trio with strongest evidence for each species conveniently
# involves all three (i.e. Pmac-Kdik-Smon). Mongo is involved in both P1-P3 and P2-P3 comparisons, but
# the strongest for mongo is between dikume and mongo (positive values). 
# Thus, we will use negative values for pungu-mongo

# When looking at P2-P3 - we will look for values > 0, when looking at P1-P3, values < 0
# Note that the strongest signal of introgression for mongo is already represented by the pmac/mongo trio
# Pmac 
# - Sympatric - P1-P3: spp1 = 'Pungu maclareni', spp3 = 'Stomatepia mongo'
pmac.sym.pre <- 'Pungu_maclareni_Konia_dikume_Stomatepia_mongo_localFstats_Dinvestigate'
# - Allopatric - P1-P3: spp1 = 'Pungu maclareni', spp3 = 'S. galileus CRM'
pmac.allo.pre <- 'Pungu_maclareni_Konia_dikume_Sarotherodon_galilaeus_CRM_localFstats_Dinvestigate'

# Kdik 
# - Sympatric - Same trio as pungu, but now between Kdik and Smon (P2-P3):
# spp1 = 'Konia dikume', spp3 = 'Stomatepia mongo'
kdik.sym.pre <- pmac.sym.pre
# - Allopatric - Same trio as pungu, but now between Kdik and SgalCRM (P2-P3): 
# spp1 = 'Konia dikume', spp3 = 'S. galileus CRM'
kdik.allo.pre <- pmac.allo.pre

# Smon Sympatric is already represented by Pmac sym - between Pmac & Smon
# - Allopatric (P2-P3): spp1 = 'Stomatepia mongo', spp3 = 'S. galileus CRM'
smon.allo.pre <- 'Stomatepia_mongo_Konia_dikume_Sarotherodon_galilaeus_CRM_localFstats_Dinvestigate'

spps <- c("Pmac.sym", "Pmac.allo", "Kdik.sym", "Kdik.allo", "Smon.allo")
spp.pairs <- list(Pmac.sym = c("Pungu maclareni", "Stomatepia mongo"), 
                  Pmac.allo = c("Pungu maclareni", "S. galilaeus CRM"),
                  Kdik.sym = c("Konia dikume", "Stomatepia mongo"), 
                  Kdik.allo = c("Konia dikume", "S. galilaeus CRM"),
                  Smon.allo = c("Stomatepia mongo", "S. galilaeus CRM"))
spp.fpaths <- c(pmac.sym.pre, pmac.allo.pre, kdik.sym.pre, kdik.allo.pre, smon.allo.pre)

chr.sizes <- read.table('../../../../chrom.sizes')
chr.sizes$V2 <- gsub(chr.sizes$V2, pattern = "LG", replacement = "")
chr.winds <- list()
for(i in 1:length(chr.sizes$V2)){
  chr.winds[[i]] <- seq(0, chr.sizes$V3[i], by = 5000)
}
names(chr.winds) <- chr.sizes$V2


source('../../../scripts/PlotSummarizeDstats.R')
windows <- downsamp.dstats(spps, spp.fpaths)

# First with pungu
pmac.sym <- 
  plt.dinvestigate(spp1 = windows$AllRes$Pmac.sym[1], 
                   spp2 = windows$AllRes$Pmac.sym[2], 
                   res = windows$DownsampRes$Pmac.sym, 
                   sim.prefix = pmac.sym.pre, thresh = 3,
                   pair = 'P1-P3', stat = 'f_dM')

pmac.allo <- 
  plt.dinvestigate(spp1 = windows$AllRes$Pmac.allo[1], 
                   spp2 = windows$AllRes$Pmac.allo[2], 
                   res = windows$DownsampRes$Pmac.allo, 
                   sim.prefix = pmac.allo.pre, thresh = 3,
                   pair = 'P1-P3', stat = 'f_dM')

# Then with dikume
kdik.sym <- 
  plt.dinvestigate(spp1 = windows$AllRes$Kdik.sym[1], 
                   spp2 = windows$AllRes$Kdik.sym[2], 
                   res = windows$DownsampRes$Kdik.sym, 
                   sim.prefix = kdik.sym.pre, thresh = 3,
                   pair = 'P2-P3', stat = 'f_dM')
kdik.allo <- 
  plt.dinvestigate(spp1 = windows$AllRes$Kdik.allo[1], 
                   spp2 = windows$AllRes$Kdik.allo[2], 
                   res = windows$DownsampRes$Kdik.allo, 
                   sim.prefix = kdik.allo.pre, thresh = 3,
                   pair = 'P2-P3', stat = 'f_dM')
# And lastly with mongo
# Note that the strongest signal of introgression for mongo is already represented by the pmac/mongo trio
smon.allo <- 
  plt.dinvestigate(spp1 = windows$AllRes$Smon.allo[1], 
                   spp2 = windows$AllRes$Smon.allo[2], 
                   res = windows$DownsampRes$Smon.allo, 
                   sim.prefix = smon.allo.pre, thresh = 3,
                   pair = 'P2-P3', stat = 'f_dM')

# Now, plot together, within species
#Pmac 
pmac.fdm <- 
  plot_grid(pmac.sym$fdm.plot + theme(axis.text.x = element_blank(), 
                                      axis.title.x = element_blank()), 
            pmac.allo$fdm.plot, ncol = 1, rel_heights = c(0.85,1))
ggsave(pmac.fdm, filename = 'Pmac-fDM.png', width = 12, height = 5.5)
ggsave(pmac.fdm, filename = 'Pmac-fDM.pdf', width = 12, height = 5.5)

pmac.df <- 
  plot_grid(pmac.sym$df.plot + theme(axis.text.x = element_blank(), 
                                     axis.title.x = element_blank()), 
            pmac.allo$df.plot, ncol = 1, rel_heights = c(0.85,1))
ggsave(pmac.df, filename = 'Pmac-dF.png', width = 12, height = 5.5)
ggsave(pmac.df, filename = 'Pmac-dF.pdf', width = 12, height = 5.5)

# Kdik
kdik.fdm <- 
  plot_grid(kdik.sym$fdm.plot + theme(axis.text.x = element_blank(), 
                                      axis.title.x = element_blank()), 
            kdik.allo$fdm.plot, ncol = 1, rel_heights = c(0.85,1))
ggsave(kdik.fdm, filename = 'Kdik-fDM.png', width = 12, height = 5.5)
ggsave(kdik.fdm, filename = 'Kdik-fDM.pdf', width = 12, height = 5.5)

kdik.df <- 
  plot_grid(kdik.sym$df.plot + theme(axis.text.x = element_blank(), 
                                     axis.title.x = element_blank()), 
            kdik.allo$df.plot, ncol = 1, rel_heights = c(0.85,1))
ggsave(kdik.df, filename = 'Kdik-dF.png', width = 12, height = 5.5)
ggsave(kdik.df, filename = 'Kdik-dF.pdf', width = 12, height = 5.5)

# Smon
ggsave(smon.allo$fdm.plot, filename = 'Smon-fDM.png', width = 12, height = 2.75)
ggsave(smon.allo$fdm.plot, filename = 'Smon-fDM.pdf', width = 12, height = 2.75)

ggsave(smon.allo$df.plot, filename = 'Smon-dF.png', width = 12, height = 2.75)
ggsave(smon.allo$df.plot, filename = 'Smon-dF.pdf', width = 12, height = 2.75)


# Put them all together. 
plot_grid(pmac.sym$fdm.plot, pmac.allo$fdm.plot,
          kdik.sym$fdm.plot, kdik.allo$fdm.plot,
          smon.allo$fdm.plot)

# Now, compile and output the significant snps for each comparison
signif.d <- pmac.smon$SiteDstats
signif.d$scaffCol <- as.character(signif.d$scaffCol)
signif.d$Pair <- "Pmac_Smon"
signif.d <- signif.d[which(signif.d$scaffCol %in% c('C', 'D')),]
signif.d$scaffCol <- 
  gsub(signif.d$scaffCol, pattern = 'C', replacement = 3) %>%
  gsub(., pattern = 'D', replacement = 4) %>%
  as.numeric(.)
signif.d <- signif.d[,c(5,1,2:4,7:9)]
colnames(signif.d) <- 
  c('Chromosome', 'Chromosome', 'WindowStart', 'WindowEnd', 
    'Fdm', 'Z-threshold', 'ObsChromZscore')


