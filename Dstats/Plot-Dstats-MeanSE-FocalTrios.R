library(tidyverse)
library(ggplot2)

# Plot the mean and SE of f_dM and d_f for each signififcant trio
# We are interested in the relative magnitude of di-statistics within 
# Barombi Mbo, as compared to among Barombi Mbi and Riverine Pops

setwd('~/Dropbox/Research/Martin-Berkeley-Postdoc/cichlids/cameroon/Onil_UMD/Dsuite/Results/100Snp_Slide25/FocalTrios/')

fpaths <- list.files(pattern="Window_100_25")

# Initialize the result dataframes
# We will create two - one for positive (i.e. P2/P3) values, and one for negative (P1/P3) values
# once done, we will combine for plotting
# Negative values (P1/P3) are not informative for Fd, so we will include in datafame, but not use.
Dstats.pos <- 
  data.frame(Trio = rep(NA, length(fpaths)), P1 = rep(NA, length(fpaths)), 
             P2 = rep(NA, length(fpaths)), P3 = rep(NA, length(fpaths)), 
             Type = rep(NA, length(fpaths)), D.mean = rep(NA, length(fpaths)), 
             D.SE = rep(NA, length(fpaths)), D.SD = rep(NA, length(fpaths)), 
             D.CI.lower = rep(NA, length(fpaths)), D.CI.upper = rep(NA, length(fpaths)), 
             f_d.mean = rep(NA, length(fpaths)), f_d.SE = rep(NA, length(fpaths)), 
             f_d.SD = rep(NA, length(fpaths)), f_d.CI.lower = rep(NA, length(fpaths)), 
             f_d.CI.upper = rep(NA, length(fpaths)), f_dM.mean = rep(NA, length(fpaths)),
             f_dM.SE = rep(NA, length(fpaths)), f_dM.SD = rep(NA, length(fpaths)), 
             f_dM.CI.lower = rep(NA, length(fpaths)), f_dM.CI.upper = rep(NA, length(fpaths)), 
             d_f.mean = rep(NA, length(fpaths)), d_f.SE = rep(NA, length(fpaths)), 
             d_f.SD = rep(NA, length(fpaths)), d_f.CI.lower = rep(NA, length(fpaths)), 
             d_f.CI.upper = rep(NA, length(fpaths)))
Dstats.neg <- Dstats.pos

# function to calculate SE
se <- function(x) sqrt(var(x)/length(x))

# A little function to help extracting species names
getSppNames <- 
  function(fpath){
    # Strip extraneous stuff from the filename to get species identies
    trio <- gsub(fpath, pattern = "_local.*", replacement = "")
    
    trio <- 
      gsub(trio, pattern = "_Konia", replacement = "-Konia") %>%
      gsub(., pattern = "_Pungu", replacement = "-Pungu") %>%
      gsub(., pattern = "_Stomatepia", replacement = "-Stomatepia") %>%
      gsub(., pattern = "_Sarotherodon", replacement = "-Sarotherodon") %>%
      gsub(., pattern = "_Myaka", replacement = "-Myaka")
    
    spp <- strsplit(trio, "-")
    
    # Clean for plotting
    trio <- 
      gsub(trio, pattern = "Stomatepia_mongo", replacement = "Smon") %>%
      gsub(., pattern = "Stomatepia_mariae", replacement = "Smar") %>%
      gsub(., pattern = "Sarotherodon_steinbachi", replacement = "Sste") %>%
      gsub(., pattern = "Sarotherodon_carolii", replacement = "Scar") %>%
      gsub(., pattern = "Stomatepia_pindu", replacement = "Spin") %>%
      gsub(., pattern = "Konia_dikume", replacement = "Kdik") %>%
      gsub(., pattern = "Konia_eisentrauti", replacement = "Keis") %>%
      gsub(., pattern = "Pungu_maclareni", replacement = "Pmac") %>%
      gsub(., pattern = "Sarotherodon_knauerae", replacement = "Pkna") %>%
      gsub(., pattern = "Sarotherodon_galilaeus_CRM", replacement = "Sgal_CRM") %>%
      gsub(., pattern = "Sarotherodon_galilaeus_MM", replacement = "Sgal_MM") %>%
      gsub(., pattern = "Myaka_myaka", replacement = "Mmya")    
    return(list(trio, spp))
  }

per.wind.res <- list()

# Will need to summarize for both positive and negative values (i.e., for P2/P3, and P1/P3 comparisons)
for(i in 1:length(fpaths)){
  print(paste("Processing Trio", i, "of", length(fpaths)))
  res <- read.table(fpaths[i], header = T)
  print(paste0("Processing a total of ", nrow(res), " windows."))
  # Get the trio of species
  trio <- getSppNames(fpaths[[i]])
  print(paste0("Trio includes P1: ", trio[[2]][[1]][1], 
               ", P2: ", trio[[2]][[1]][2], 
               ", and P3: ", trio[[2]][[1]][3]))
  
  # Put into the species identities into the dataframe
  # First postive (P2/P3)
  Dstats.pos[i,'Trio'] <- trio[[1]][[1]]
  Dstats.pos[i,2:4] <- trio[[2]][[1]]
  if(trio[[2]][[1]][3] %in% c('Sarotherodon_knauerae',
                              'Sarotherodon_galilaeus_MM',
                              'Sarotherodon_galilaeus_CRM')){
    Dstats.pos[i,'Type'] <- "Allopatric"
  }else{
    Dstats.pos[i,'Type'] <- "Sympatric"
  }
  
  # Sometimes f_d appears to be estimated incorrectly, with values > 1 or < 0. 
  # Check if these values exist, and remove
  too.big <- which(res$f_d > 1 | res$f_d < -1)
  if(length(too.big > 0)){
    res <- res[-too.big,]
  }
  print(paste0("Removed ", length(too.big), " windows with values of f_d that were too large ( -1 > f_d > 1"))
  
  Dstats.pos[i,'D.mean'] <- mean(res$D[which(res$D > 0)])
  Dstats.pos[i,'D.SE'] <- se(res$D[which(res$D > 0)])
  Dstats.pos[i,'D.SD'] <- sd(res$D[which(res$D > 0)])
  Dstats.pos[i,'D.CI.lower'] <- quantile(res$D[which(res$D > 0)], probs = 0.025)
  Dstats.pos[i,'D.CI.upper'] <- quantile(res$D[which(res$D > 0)], probs = 0.925)
  
  Dstats.pos[i,'f_d.mean'] <- mean(res$f_d[which(res$f_d > 0)])
  Dstats.pos[i,'f_d.SE'] <- se(res$f_d[which(res$f_d > 0)])
  Dstats.pos[i,'f_d.SD'] <- sd(res$f_d[which(res$f_d > 0)])
  Dstats.pos[i,'f_d.CI.lower'] <- quantile(res$f_d[which(res$f_d > 0)], probs = 0.025)
  Dstats.pos[i,'f_d.CI.upper'] <- quantile(res$f_d[which(res$f_d > 0)], probs = 0.925)
  
  Dstats.pos[i,'f_dM.mean'] <- mean(res$f_dM[which(res$f_dM > 0)])
  Dstats.pos[i,'f_dM.SE'] <- se(res$f_dM[which(res$f_dM > 0)])
  Dstats.pos[i,'f_dM.SD'] <- sd(res$f_dM[which(res$f_dM > 0)])
  Dstats.pos[i,'f_dM.CI.lower'] <- quantile(res$f_dM[which(res$f_dM > 0)], probs = 0.025)
  Dstats.pos[i,'f_dM.CI.upper'] <- quantile(res$f_dM[which(res$f_dM > 0)], probs = 0.925)
  
  Dstats.pos[i,'d_f.mean'] <- mean(res$d_f[which(res$d_f > 0)])
  Dstats.pos[i,'d_f.SE'] <- se(res$d_f[which(res$d_f > 0)])
  Dstats.pos[i,'d_f.SD'] <- sd(res$d_f[which(res$d_f > 0)])
  Dstats.pos[i,'d_f.CI.lower'] <- quantile(res$d_f[which(res$d_f > 0)], probs = 0.025)
  Dstats.pos[i,'d_f.CI.upper'] <- quantile(res$d_f[which(res$d_f > 0)], probs = 0.925)
  
  # Then negative (P1/P3)
  Dstats.neg[i,'Trio'] <- trio[[1]][[1]]
  Dstats.neg[i,2:4] <- trio[[2]][[1]]
  if(trio[[2]][[1]][3] %in% c('Sarotherodon_knauerae',
                              'Sarotherodon_galilaeus_MM', 
                              'Sarotherodon_galilaeus_CRM')){
    Dstats.neg[i,'Type'] <- "Allopatric"
  }else{
    Dstats.neg[i,'Type'] <- "Sympatric"
  }
  Dstats.neg[i,'D.mean'] <- mean(res$D[which(res$D < 0)])
  Dstats.neg[i,'D.SE'] <- se(res$D[which(res$D < 0)])
  Dstats.neg[i,'D.SD'] <- sd(res$D[which(res$D < 0)])
  Dstats.neg[i,'D.CI.lower'] <- quantile(res$D[which(res$D < 0)], probs = 0.025)
  Dstats.neg[i,'D.CI.upper'] <- quantile(res$D[which(res$D < 0)], probs = 0.925)
  
  Dstats.neg[i,'f_d.mean'] <- mean(res$f_d[which(res$f_d < 0)])
  Dstats.neg[i,'f_d.SE'] <- se(res$f_d[which(res$f_d < 0)])
  Dstats.neg[i,'f_d.SD'] <- sd(res$f_d[which(res$f_d < 0)])
  Dstats.neg[i,'f_d.CI.lower'] <- quantile(res$f_d[which(res$f_d < 0)], probs = 0.025)
  Dstats.neg[i,'f_d.CI.upper'] <- quantile(res$f_d[which(res$f_d < 0)], probs = 0.925)
  
  Dstats.neg[i,'f_dM.mean'] <- mean(res$f_dM[which(res$f_dM < 0)])
  Dstats.neg[i,'f_dM.SE'] <- se(res$f_dM[which(res$f_dM < 0)])
  Dstats.neg[i,'f_dM.SD'] <- sd(res$f_dM[which(res$f_dM < 0)])
  Dstats.neg[i,'f_dM.CI.lower'] <- quantile(res$f_dM[which(res$f_dM < 0)], probs = 0.025)
  Dstats.neg[i,'f_dM.CI.upper'] <- quantile(res$f_dM[which(res$f_dM < 0)], probs = 0.925)
  
  Dstats.neg[i,'d_f.mean'] <- mean(res$d_f[which(res$d_f < 0)])
  Dstats.neg[i,'d_f.SE'] <- se(res$d_f[which(res$d_f < 0)])
  Dstats.neg[i,'d_f.SD'] <- sd(res$d_f[which(res$d_f < 0)])
  Dstats.neg[i,'d_f.CI.lower'] <- quantile(res$d_f[which(res$d_f < 0)], probs = 0.025)
  Dstats.neg[i,'d_f.CI.upper'] <- quantile(res$d_f[which(res$d_f < 0)], probs = 0.925)
  
  # Now add the per-window d-stats
  per.wind.res[[i]] <- res
  names(per.wind.res)[i] <- trio[[1]]
}

# Now combine
Dstats.pos$Pair <- "P2-P3"
Dstats.neg$Pair <- "P1-P3"
Dstats <- rbind(Dstats.pos, Dstats.neg)

write.table(Dstats, file = 'D-Stats_Mean-SE_SignifTrios_100Snp_Slide25.tsv', 
            sep ="\t", row.names = T, quote = F)
Dstats <- read.table('D-Stats_Mean-SE_SignifTrios_100Snp_Slide25.tsv')

# Some stuff to make it behave when plotting
Dstats$Type <- 
  factor(Dstats$Type, levels = c("Sympatric", "Allopatric"))

Dstats$Trio <- 
  factor(Dstats$Trio, levels = Dstats$Trio[order(Dstats$d_f.mean[which(Dstats$d_f.mean > 0)])])

# order them neatly
Dstats$Trio <- 
  factor(Dstats$Trio, 
         levels = c("Pmac-Smon-Sgal_CRM", "Kdik-Smon-Sgal_CRM",
                    "Smon-Kdik-Sgal_CRM", "Pmac-Kdik-Sgal_CRM",
                    "Smon-Pmac-Sgal_CRM", "Kdik-Pmac-Sgal_CRM",
                    "Smon-Pmac-Kdik", "Smon-Kdik-Pmac",
                    "Pmac-Kdik-Smon", "Pmac-Smon-Kdik",
                    "Kdik-Smon-Pmac","Kdik-Pmac-Smon"))

WindStats <- 
  data.frame(
    "Trio" = c(rep(names(per.wind.res)[1], length = nrow(per.wind.res[[1]])),
               rep(names(per.wind.res)[2], length = nrow(per.wind.res[[2]])),
               rep(names(per.wind.res)[3], length = nrow(per.wind.res[[3]])),
               rep(names(per.wind.res)[4], length = nrow(per.wind.res[[4]])),
               rep(names(per.wind.res)[5], length = nrow(per.wind.res[[5]])),
               rep(names(per.wind.res)[6], length = nrow(per.wind.res[[6]])),
               rep(names(per.wind.res)[7], length = nrow(per.wind.res[[7]])),
               rep(names(per.wind.res)[8], length = nrow(per.wind.res[[8]])),
               rep(names(per.wind.res)[9], length = nrow(per.wind.res[[9]])),
               rep(names(per.wind.res)[10], length = nrow(per.wind.res[[10]])),
               rep(names(per.wind.res)[11], length = nrow(per.wind.res[[11]])),
               rep(names(per.wind.res)[12], length = nrow(per.wind.res[[12]]))),
    "D" = c(per.wind.res[[1]]$D, per.wind.res[[2]]$D, per.wind.res[[3]]$D,
            per.wind.res[[4]]$D, per.wind.res[[5]]$D, per.wind.res[[6]]$D,
            per.wind.res[[7]]$D, per.wind.res[[8]]$D, per.wind.res[[9]]$D,
            per.wind.res[[10]]$D, per.wind.res[[11]]$D, per.wind.res[[12]]$D),
    "f_d" = c(per.wind.res[[1]]$f_d, per.wind.res[[2]]$f_d, per.wind.res[[3]]$f_d,
            per.wind.res[[4]]$f_d, per.wind.res[[5]]$f_d, per.wind.res[[6]]$f_d,
            per.wind.res[[7]]$f_d, per.wind.res[[8]]$f_d, per.wind.res[[9]]$f_d,
            per.wind.res[[10]]$f_d, per.wind.res[[11]]$f_d, per.wind.res[[12]]$f_d),
    "f_dM" = c(per.wind.res[[1]]$f_dM, per.wind.res[[2]]$f_dM, per.wind.res[[3]]$f_dM,
            per.wind.res[[4]]$f_dM, per.wind.res[[5]]$f_dM, per.wind.res[[6]]$f_dM,
            per.wind.res[[7]]$f_dM, per.wind.res[[8]]$f_dM, per.wind.res[[9]]$f_dM,
            per.wind.res[[10]]$f_dM, per.wind.res[[11]]$f_dM, per.wind.res[[12]]$f_dM),
    "d_f" = c(per.wind.res[[1]]$d_f, per.wind.res[[2]]$d_f, per.wind.res[[3]]$d_f,
            per.wind.res[[4]]$d_f, per.wind.res[[5]]$d_f, per.wind.res[[6]]$d_f,
            per.wind.res[[7]]$d_f, per.wind.res[[8]]$d_f, per.wind.res[[9]]$d_f,
            per.wind.res[[10]]$d_f, per.wind.res[[11]]$d_f, per.wind.res[[12]]$d_f)
  )

WindStats$Type <- NA

allo <- c("Pmac-Smon-Sgal_CRM", "Kdik-Smon-Sgal_CRM",
          "Smon-Kdik-Sgal_CRM", "Pmac-Kdik-Sgal_CRM",
          "Smon-Pmac-Sgal_CRM", "Kdik-Pmac-Sgal_CRM")
symp <- c("Smon-Pmac-Kdik", "Smon-Kdik-Pmac",
          "Pmac-Kdik-Smon", "Pmac-Smon-Kdik",
          "Kdik-Smon-Pmac","Kdik-Pmac-Smon")
WindStats$Type[which(WindStats$Trio %in% allo)] <- "Allopatric"
WindStats$Type[which(WindStats$Trio %in% symp)] <- "Sympatric"

WindStats$Trio <- 
  factor(WindStats$Trio, 
         levels = c("Pmac-Smon-Sgal_CRM", "Kdik-Smon-Sgal_CRM",
                    "Smon-Kdik-Sgal_CRM", "Pmac-Kdik-Sgal_CRM",
                    "Smon-Pmac-Sgal_CRM", "Kdik-Pmac-Sgal_CRM",
                    "Smon-Pmac-Kdik", "Smon-Kdik-Pmac",
                    "Pmac-Kdik-Smon", "Pmac-Smon-Kdik",
                    "Kdik-Smon-Pmac","Kdik-Pmac-Smon"))
WindStats$Type <- 
  factor(WindStats$Type, levels = c("Sympatric", "Allopatric"))
WindStats$D.Pair <- NA
WindStats$D.Pair[which(WindStats$D > 0)] <- "P2-P3"
WindStats$D.Pair[which(WindStats$D < 0)] <- "P1-P3"
WindStats$f_d.Pair <- NA
WindStats$f_d.Pair[which(WindStats$f_d > 0)] <- "P2-P3"
WindStats$f_d.Pair[which(WindStats$f_d < 0)] <- "P1-P3"
WindStats$f_dM.Pair <- NA
WindStats$f_dM.Pair[which(WindStats$f_dM > 0)] <- "P2-P3"
WindStats$f_dM.Pair[which(WindStats$f_dM < 0)] <- "P1-P3"
WindStats$d_f.Pair <- NA
WindStats$d_f.Pair[which(WindStats$d_f > 0)] <- "P2-P3"
WindStats$d_f.Pair[which(WindStats$d_f < 0)] <- "P1-P3"
# Are sympatric trios characterized by greater amounts of introgression?
#p <- t.test(Dstats$d_f.mean[which(Dstats$Type == "Sympatric")], Dstats$d_f.mean[which(Dstats$Type == "Allopatric")])$p.value
df.plot <- 
  ggplot(data = Dstats, aes(y = Trio, x = d_f.mean)) +
  geom_violin(data = WindStats, aes(y = Trio, x = d_f, fill = Type), 
              alpha = 0.5, draw_quantiles = c(0.25, 0.75)) + 
  scale_fill_manual(values = c('#74add1', '#f46d43')) + 
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(xmin = d_f.mean-(2*d_f.SE), xmax = d_f.mean+(2*d_f.SE)), size = 1, color = 'black') +
  xlab(expression(italic('df'))) + 
  geom_point(size = 4, color = 'black') + 
  theme_bw(base_size = 20) + 
  coord_cartesian(xlim = c(-0.3, 0.3)) +
  theme(legend.position = 'top', 
        legend.title = element_blank())
df.plot
ggsave(df.plot, filename = "d_f-GenomeWide-MeanSE_100Snp_Slide25-Focal.pdf", height = 18, width = 12)

Dstats$Trio <- 
  factor(Dstats$Trio, levels = Dstats$Trio[order(Dstats$f_dM.mean[which(Dstats$f_dM.mean > 0)])])
fdm.plot <- 
  ggplot(data = Dstats, aes(y = Trio, x = f_dM.mean)) +
  geom_violin(data = WindStats, aes(y = Trio, x = f_dM, fill = Type), 
              alpha = 0.5, draw_quantiles = c(0.25, 0.75)) + 
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(xmin = f_dM.mean-(2*f_dM.SE), xmax = f_dM.mean+(2*f_dM.SE)), size = 1.25) +
  scale_fill_manual(values = c('#74add1', '#f46d43')) + 
  xlab(expression(italic(f)[dM])) + 
  geom_point(size = 4, color = 'black') + 
  # annotate('text', x = 0, y = 50, hjust = -0.1, size = 8,
  #          label = paste0("t-test = ", round(p, 3))) +
  theme_bw(base_size = 20) + 
  coord_cartesian(xlim = c(-0.2, 0.2)) +
  theme(legend.position = 'top', 
        legend.title = element_blank())
ggsave(fdm.plot, filename = "fdM-GenomeWide-MeanSE_100Snp_Slide25-Focal.pdf", height = 18, width = 12)

# Dstats$Trio <- 
#   factor(Dstats$Trio, levels = Dstats$Trio[order(Dstats$D.mean[which(Dstats$D.mean > 0)])])
# p <- t.test(Dstats$D.mean[which(Dstats$Type == "Sympatric")], Dstats$D.mean[which(Dstats$Type == "Allopatric")])$p.value
d.plot <- 
  ggplot(data = Dstats, aes(y = Trio, x = D.mean)) +
  geom_violin(data = WindStats, aes(y = Trio, x = D, fill = Type), 
              alpha = 0.5, draw_quantiles = c(0.25, 0.75)) + 
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(xmin = D.mean-(2*D.SE), xmax = D.mean+(2*D.SE)), size = 1.25) +
  scale_fill_manual(values = c('#74add1', '#f46d43')) + 
  xlab(expression(italic('D'))) + 
  geom_point(size = 4, color = 'black') + 
  # annotate('text', x = 0, y = 50, hjust = -0.1, size = 8,
  #          label = paste0("t-test = ", round(p, 3))) +
  theme_bw(base_size = 20) + 
  coord_cartesian(xlim = c(-0.4, 0.4)) +
  theme(legend.position = 'top', 
        legend.title = element_blank())
ggsave(d.plot, filename = "D-GenomeWide-MeanSE_100Snp_Slide25-Focal.pdf", height = 18, width = 12)

# p <- t.test(Dstats$f_d.mean[which(Dstats$Type == "Sympatric")], Dstats$f_d.mean[which(Dstats$Type == "Allopatric")])$p.value
# fdStats$Trio <- 
#   factor(fdStats$Trio, levels = fdStats$Trio[order(fdStats$f_d.mean)])

fd.plot <- 
  ggplot(data = Dstats, aes(y = Trio, x = f_d.mean)) +
  geom_violin(data = WindStats, aes(y = Trio, x = f_d, fill = Type), 
              alpha = 0.5, draw_quantiles = c(0.25, 0.75)) + 
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(xmin = f_d.mean-(2*f_d.SE), xmax = f_d.mean+(2*f_d.SE)), size = 1.25) +
  scale_fill_manual(values = c('#74add1', '#f46d43')) + 
  xlab(expression(italic(f)[d])) + 
  geom_point(size = 4, color = 'black') + 
  # annotate('text', x = 0, y = 50, hjust = -0.1, size = 8,
  #          label = paste0("t-test = ", round(p, 3))) +
  theme_bw(base_size = 20) + 
  coord_cartesian(xlim = c(-0.3, 0.3)) +
  theme(legend.position = 'top', 
        legend.title = element_blank())
ggsave(fd.plot, filename = "fd-GenomeWide-MeanSE_100Snp_Slide25-Focal.pdf", height = 18, width = 12)

#################################################################################
# Plot the subset for which only the focal species are in P2 or P3
focal.spp <- c('Konia_dikume', 'Pungu_maclareni', 'Stomatepia_mongo')
Dstats <- Dstats[which(Dstats$P1 %in%focal.spp | Dstats$P2 %in% focal.spp),]
fdStats <- fdStats[which(fdStats$P1 %in%focal.spp | fdStats$P2 %in% focal.spp),]

# Some stuff to make it behave when plotting
Dstats$Type <- 
  factor(Dstats$Type, levels = c("Sympatric", "Allopatric"))
Dstats$Trio <- 
  factor(Dstats$Trio, levels = Dstats$Trio[order(Dstats$d_f.mean[which(Dstats$d_f.mean > 0)])])

df.plot <- 
  ggplot(data = Dstats, aes(y = Trio, x = d_f.mean, color = Type)) +
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(xmin = d_f.mean-(2*d_f.SE), xmax = d_f.mean+(2*d_f.SE)), size = 1.25) +
  scale_color_manual(values = c('#74add1', '#f46d43')) +
  xlab(expression(italic('df'))) + 
  geom_point(size = 4) + 
  theme_bw(base_size = 20) + 
  theme(legend.position = 'top', 
        legend.title = element_blank())
df.plot
ggsave(df.plot, filename = "d_f-GenomeWide-MeanSE-FocalSpp_100Snp_Slide25-Focal.pdf", height = 10, width = 10)

Dstats$Trio <- 
  factor(Dstats$Trio, levels = Dstats$Trio[order(Dstats$f_dM.mean[which(Dstats$f_dM.mean > 0)])])
fdm.plot <- 
  ggplot(data = Dstats, aes(y = Trio, x = f_dM.mean, color = Type)) +
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(xmin = f_dM.mean-(2*f_dM.SE), xmax = f_dM.mean+(2*f_dM.SE)), size = 1.25) +
  scale_color_manual(values = c('#74add1', '#f46d43')) +
  xlab(expression(italic(f)[dM])) + 
  geom_point(size = 4) + 
  theme_bw(base_size = 20) + 
  theme(legend.position = 'top', 
        legend.title = element_blank())
ggsave(fdm.plot, filename = "fdM-GenomeWide-MeanSE-FocalSpp_100Snp_Slide25-Focal.pdf", height = 10, width = 10)

Dstats$Trio <- 
  factor(Dstats$Trio, levels = Dstats$Trio[order(Dstats$D.mean[which(Dstats$D.mean > 0)])])
d.plot <- 
  ggplot(data = Dstats, aes(y = Trio, x = D.mean, color = Type)) +
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(xmin = D.mean-(2*D.SE), xmax = D.mean+(2*D.SE)), size = 1.25) +
  scale_color_manual(values = c('#74add1', '#f46d43')) +
  xlab(expression(italic('D'))) + 
  geom_point(size = 4) + 
  theme_bw(base_size = 20) + 
  theme(legend.position = 'top', 
        legend.title = element_blank())
ggsave(d.plot, filename = "D-GenomeWide-MeanSE-FocalSpp_100Snp_Slide25-Focal.pdf", height = 10, width = 10)

fd.plot <- 
  ggplot(data = fdStats, aes(y = Trio, x = f_d.mean, color = Type)) +
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(xmin = f_d.mean-(2*f_d.SE), xmax = f_d.mean+(2*f_d.SE)), size = 1.25) +
  scale_color_manual(values = c('#74add1', '#f46d43')) +
  xlab(expression(italic(f)[d])) + 
  geom_point(size = 4) + 
  theme_bw(base_size = 20) + 
  theme(legend.position = 'top', 
        legend.title = element_blank())
ggsave(fd.plot, filename = "fd-GenomeWide-MeanSE-FocalSpp_100Snp_Slide25-Focal.pdf", height = 10, width = 10)



# Now, plot Fdm and Df together
Dstats$Trio <-  
  factor(Dstats$Trio, levels = Dstats$Trio[order(Dstats$f_dM.mean[which(Dstats$f_dM.mean > 0)])])

df.plot <- 
  ggplot(data = Dstats, aes(y = Trio, x = d_f.mean, color = Type)) +
  geom_vline(xintercept = 0, lty = 2, alpha = 0.5, size = 0.5) +
  geom_errorbar(aes(xmin = d_f.mean-(2*d_f.SE), xmax = d_f.mean+(2*d_f.SE)), size = 1.25) +
  scale_color_manual(values = c('#74add1', '#f46d43')) +
  xlab(expression(italic('df'))) + 
  geom_point(size = 4) + 
  theme_bw(base_size = 20) + 
  theme(legend.position = 'top', 
        legend.title = element_blank())
df.plot

focal_fdm_df.p <- 
  plot_grid(fdm.plot, df.plot + 
              theme(axis.text.y = element_blank(), 
                    axis.ticks.y = element_blank()),
            rel_widths = c(1, 0.725))

ggsave(focal_fdm_df.p, filename = "fdm-df-GenomeWide-MeanSE-FocalSpp_100Snp_Slide25-Focal.pdf", height = 10, width = 15)



