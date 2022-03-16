# downsample such that we only include one window per 5kb
downsamp.dstats <- 
  function(spps, spp.fpaths){
    spp.res <- list()
    spp.res.downsamp <- list()
    keeper.winds <- c()
    
    for(spp in 1:length(spps)){
      Species <- spps[spp]
      print(paste0("Beginning ", Species))
      
      all.chr <- read.table(list.files(pattern = spp.fpaths[spp]), header = T)
      tmp$Chromosome <- NA
      
      print(paste0("There are a total of ", nrow(all.chr), " windows."))
      
      n = 0
      for(i in 1:nrow(chr.sizes)){
        print(paste0("Working on chromosome ", chr.sizes$V2[i]))
        
        n <- n+1
        tmp <- all.chr[which(all.chr$chr == chr.sizes$V1[i]),]
        tmp$Chromosome <- chr.sizes$V2[i]
        tmp$ID <- paste0(tmp$Chromosome, "_", 1:nrow(tmp))
        
        
        tmp$windowMid <- round(rowMeans(tmp[,2:3]))
        if(n == 1){
          tmp$cumPos <- tmp$windowMid
          tmp$cumStart <- tmp$windowStart
          tmp$cumEnd <- tmp$windowEnd
          res <- tmp
        }else{
          tmp$cumPos <- tmp$windowMid + max(res$cumPos)
          tmp$cumStart <- tmp$windowStart + max(res$cumStart)
          tmp$cumEnd <- tmp$windowEnd + max(res$cumEnd)
          res <- rbind(res, tmp)
        }
        
        winds.keep <- 
          data.frame(start = chr.winds[[i]]+1,
                     end = chr.winds[[i]]+5000, 
                     NumWinds = NA,
                     KeepWind = NA,
                     ID = paste0(paste0(chr.sizes$V2[i], "_", 
                                        1:length(chr.winds[[i]]))))
        
        for(w in 1:nrow(winds.keep)){
          focal <- which(tmp$windowMid >= winds.keep$start[w] & tmp$windowMid <= winds.keep$end[w])
          if(length(focal) == 0){
            winds.keep$NumWinds[w] <- 0
            winds.keep$KeepWind[w] <- NA
          }else{
            # This keeps the max statistic in each 5kb window. 
            #keeper <- focal[which(tmp[focal,7] == max(tmp[focal,7]))]
            # What if we just randomly sample per 5kb window. Less biased?
            # Potentially also less power. 
            # Remove any window that's already been sampled
            dups <- which(focal %in% unique(winds.keep$KeepWind))
            if(length(dups) > 0){
              focal <- focal[-dups]
            }
            keeper <- focal[sample(1:length(focal), size = 1)]
            winds.keep$NumWinds[w] <- length(focal)
            winds.keep$KeepWind[w] <- keeper
          }
        }
        
        winds.keep <- na.omit(winds.keep)
        tmp <- tmp[winds.keep$KeepWind,]
        tmp$WindStart.5k <- winds.keep$start
        tmp$WindEnd.5k <- winds.keep$end
        
        if(n == 1){
          res.downsamp <- tmp
        }else{
          res.downsamp <- rbind(res.downsamp, tmp)
        }
      }
      print(paste0("Sampling one window per 5kb led to a reduction from ", nrow(all.chr), " to ", nrow(res.downsamp), " windows."))
      
      spp.res[[spp]] <- res
      names(spp.res)[spp] <- Species
      spp.res.downsamp[[spp]] <- res.downsamp
      names(spp.res.downsamp)[spp] <- Species
    }
    results <- 
      list(
        AllRes = spp.res, 
        DownsampRes = spp.res.downsamp)
    return(results)
  }

zscore <-
  function(obs.vals, ref.vals){
    m <- mean(ref.vals)
    stdev <- sd(ref.vals)
    z <- (obs.vals - m)/stdev
    return(z)
  }

merge.peaks <- function(obs, signif.metric = c('zscore', 'p'), thresh, fstat = c('fDM', 'dF')){
  
  # This function will pull out the merged 'peaks' for the raised sweep statistics
  # We will pull out both the window positions, and calculate the mean of each 
  # summary statistic for each merged peak, for each statistic, and perform 
  # a test of how unusual that value is as compared to what was observed in 
  # coalescent simulations.
  
  for(stat in fstat){
    if(stat == "fDM"){
      statistic <- "fDM"
      color.col <- 14
      if(signif.metric == 'zscore'){
        thresh.col <- 16
      }else{
        thresh.col <- 18
      }
    }else{
      statistic = "dF"
      color.col <- 15
      if(signif.metric == 'zscore'){
        thresh.col <- 17
      }else{
        thresh.col <- 19
      }          
    }
    
    print(paste0('Working on ', statistic))
    
    if(signif.metric == 'p'){
      thresh <- 0.005
      peaks <- which(obs[,thresh.col] <= thresh)
    }else{
      thresh <- thresh
      thresh.1 <- thresh
      thresh.2 <- thresh+1
      
      peaks <- which(obs[,thresh.col] >= thresh)
    }
    
    if(nrow(na.omit(obs[peaks,])) == 0){
      hdrs <- obs[peaks,]
    }else{
      hdrs <- na.omit(obs[peaks,])
    }
    
    if(nrow(hdrs) > 0){
      hdrs$Chromosome <- as.numeric(hdrs$Chromosome)
      print(paste0(length(peaks), ' significant windows identified at thresh = ', thresh.1))
      
      # put a buffer of 25kb around each peak to account for LD decay
      merged.winds <- 
        data.frame(raw.mid = obs[peaks,]$windowMid,
                   raw.start = obs[peaks,]$WindowStart,
                   raw.end = obs[peaks,]$WindowEnd,
                   midpoint = obs[peaks,]$cumPos,
                   lower = obs[peaks,]$cumStart-25000,
                   upper = obs[peaks,]$cumEnd+25000,
                   chrom = obs[peaks,]$Chromosome,
                   merged.wind.ID = 0)
      
      for(i in 1:(nrow(merged.winds))){
        if(i == 1){
          merged.w.prev <- F
        }else{
          merged.w.prev <- merged.winds$lower[i] <= merged.winds$upper[i-1]
        }
        if(i == nrow(merged.winds)){
          merge.w.next <- F
        }else{
          merge.w.next <- merged.winds$upper[i] >= merged.winds$lower[i+1]
        }
        if(merged.w.prev == F){
          merged.winds$merged.wind.ID[i] <- 
            max(merged.winds$merged.wind.ID) + 1
        }else{
          if(merged.w.prev == T & merge.w.next == T){
            merged.winds$merged.wind.ID[i] <- 
              max(merged.winds$merged.wind.ID)
          }else{
            merged.winds$merged.wind.ID[i] <- 
              max(merged.winds$merged.wind.ID)
          }
        }
      }
      
      comb.hdrs <- 
        data.frame(Chromosome = NA, Position = NA, start = NA, end = NA, 
                   NumWinds = NA, CumPos = NA, CumulativeStart = NA, CumulativeEnd = NA, 
                   fDM = NA, dF = NA, zfDM = NA, zDF = NA, 
                   fDM.fdr = NA, dF.fdr = NA, stat = statistic)
      
      # Now fill in for all the merged peaks
      for(i in 1:max(merged.winds$merged.wind.ID)){
        focal <- merged.winds[which(merged.winds$merged.wind.ID == i),]
        comb.hdrs$Chromosome[i] <- as.numeric(unique(focal$chrom))
        comb.hdrs$Position[i] <- mean(c(min(focal$raw.start), max(focal$raw.end)))
        comb.hdrs$start[i] <- min(focal$raw.start)
        comb.hdrs$end[i] <- max(focal$raw.end)
        comb.hdrs$NumWinds[i] <- nrow(focal)
        comb.hdrs$CumPos[i] <- mean(c(min(focal$lower), max(focal$upper)))
        comb.hdrs$CumulativeStart[i] <- min(focal$lower)
        comb.hdrs$CumulativeEnd[i] <- max(focal$upper)
        
        # Pull out the stats
        focal <- obs[which(obs$cumStart >= min(focal$lower) & 
                             obs$cumEnd <= max(focal$upper)),]
        comb.hdrs$fDM[i] <- mean(focal$fDM[which(focal$fDM > 0)])
        comb.hdrs$dF[i] <- mean(focal$dF[which(focal$dF > 0)])
        
        if(i < max(merged.winds$merged.wind.ID)){
          comb.hdrs[i+1,] <- NA
          comb.hdrs$stat[i+1] <- statistic
        }
      }
      
      
      print(paste0('Allowing for a linkage buffer of 50kb (25kb flanking), this is a total of ', sum(comb.hdrs$NumWinds), ' unmerged windows'))
      print(paste0('These have been merged into ', nrow(comb.hdrs), ' unique, merged windows'))
      
      # Identify the false discovery rate for the means of each summary statistic,
      # sampling 50,000 merged windows of sizes equal to those in our observed merged 
      # window set.
      sizes <- sample(comb.hdrs$NumWinds, replace = T, 50000)
      fdm.resamp <- c()
      df.resamp <- c()
      
      print(paste0('Sample 50000 randomly merged windows of equivalent size to test whether other statistics in these peaks are exceptional.'))
      for(i in 1:length(sizes)){
        size <- sizes[i]
        
        good = F
        while(good == F){
          samp.i <- sample(1:nrow(obs), 1)
          if(samp.i+size > nrow(obs)){
            # push the window back to fit
            samp.i <- samp.i - (size - (nrow(obs) - samp.i))
          }
          samp <- obs[samp.i:(samp.i+size),]
          same.chr <- length(unique(samp$Chromosome)) == 1
          neighboring <- max(diff(samp$cumPos)) <= 25000
          if(same.chr == T & neighboring == T){good = T}
        }
        fdm.resamp <- append(fdm.resamp, mean(samp$fDM))
        df.resamp <- append(df.resamp, mean(samp$dF))
      }
      
      comb.hdrs$zfDM <- zscore(comb.hdrs$fDM, fdm.resamp)
      comb.hdrs$zDF <- zscore(comb.hdrs$dF, df.resamp)
      
      comb.hdrs$fDM.fdr <- 
        p.adjust(pnorm(comb.hdrs$zfDM, 0, 1, lower.tail = F), method = 'fdr')
      comb.hdrs$dF.fdr <- 
        p.adjust(pnorm(comb.hdrs$zDF, 0, 1, lower.tail = F), method = 'fdr')
      
      if(stat == "fDM"){
        fdm.hdrs <- comb.hdrs
      }else{
        df.hdrs <- comb.hdrs
      }
    }else{
      # No windows with summary statistics that are exceptional as compared to what we
      # observe in our simulated data.
      # Fill in with NAs to make sure it doesn't just break.
      comb.hdrs <- 
        data.frame(Chromosome = NA, Position = NA, start = NA, end = NA, 
                   NumWinds = NA, CumPos = NA, CumulativeStart = NA, CumulativeEnd = NA, 
                   fDM = NA, dF = NA, zfDM = NA, zDF = NA, 
                   fDM.fdr = NA, dF.fdr = NA, stat = statistic)
      
      if(stat == 'fDM'){
        fdm.hdrs <- comb.hdrs
      }else{
        df.hdrs <- comb.hdrs
      }
    }
  }
  if(fstat == 'fDM'){
    all.hdrs <- fdm.hdrs
  }else{
    if(fstat == 'dF'){
      all.hdrs <- df.hdrs
    }else{
      all.hdrs <- rbind(fdm.hdrs, df.hdrs)
    }
  }
  return(merged.winds = all.hdrs)
}


plt.dinvestigate <- 
  function(spp1, spp2, sim.prefix, res, pair = c('P1-P3', 'P2-P3'), stat = 'f_dM', thresh){
    # You can specify different statistics if you want - we're using d_f here
    lgs <- read.table('~/Dropbox/Research/Martin-Berkeley-Postdoc/cichlids/cameroon/Onil_UMD/Onil_UMD-Scaffold-IDs.txt')
    lgs <- lgs[c(1:2,22, 3:21),]
    res$fdm.col <- NA
    res$df.col <- NA
    res$zFDM <- NA
    res$zDF <- NA
    res$zFDM.fdr <- NA
    res$zDF.fdr <- NA
    
    res <- res[,c(1:13,16:21,14:15)]
    
    sim <- read.table(paste0('../../SimResults/100Snp_Slide25/', sim.prefix, '-Sim-Window_100_25.txt'), 
                      header = T)
    
    if(spp2 %in% c("Stomatepia mongo", "Konia dikume", "Pungu maclareni")){
      patry <- 'Barombi Mbo'
    }else{
      patry <- 'S. galilaeus CRM'
    }
    
    res$Chromosome <- 
      factor(res$Chromosome, levels = c(1:20, 22:23))
    
    a <- c(1,3,5,7,9,11,13,15,17,19,22)
    b <- c(2,4,6,8,10,12,14,16,18,20,23)
    res[which(res$Chromosome %in% a),14:15] <- 'A'
    res[which(res$Chromosome %in% b),14:15] <- 'B'
    
    # And then remove negative values if looking at  P2-P3 (or the opposite if looking at P1-P3)
    if(pair == "P1-P3"){
      res <- res[which(res[,stat] < 0),]
      res[,stat] <- abs(res[,stat])
      sim <- sim[which(sim[,stat] < 0),]
      sim[,stat] <- abs(sim[,stat])
      # alt.res <- alt.res[which(alt.res[,stat] < 0),]
      # alt.res[,stat] <- abs(alt.res[,stat])
    }else{
      res <- res[which(res[,stat] > 0),]
      sim <- sim[which(sim[,stat] > 0),]
      # alt.res <- alt.res[which(alt.res[,stat] > 0),]
    }
    
    # First without it
    # Calculate z-scores 
    # Note that we also want to compare across all comparisons - we want allopatric and 
    # sympatric comparisons to be 'treated equally'
    # IF SIMULATIONS DON'T SEEM LIKE THEY'LL WORK, COULD POSSIBLY JUST RESAMPLE FROM FDM DISTRIBUTIONS
    # TO GET A SENSE OF WHAT THE EXPECTED/CORE DISTRIBUTION IS ACROSS ALL POSSIBLE COMPARISONS?
    # all.fdm <- c(spp.res.downsamp[[1]]$f_dM, 
    #              spp.res.downsamp[[2]]$f_dM, 
    #              spp.res.downsamp[[3]]$f_dM, 
    #              spp.res.downsamp[[4]]$f_dM)
    
    res$zFDM <- zscore(obs.vals = res$f_dM, 
                       ref.vals = sim$f_dM)
    res$zDF <- zscore(obs.vals = res$d_f, 
                      ref.vals = sim$d_f)
    
    res$zFDM.fdr <- p.adjust(pnorm(res$zFDM, 0, 1, lower.tail = F), method = 'fdr')
    res$zDF.fdr <- p.adjust(pnorm(res$zDF, 0, 1, lower.tail = F), method = 'fdr')
    
    thresh.1 <- thresh
    thresh.2 <- thresh+1
    
    res$fdm.col[which(res$zFDM >= thresh.1)] <- 'C'
    res$fdm.col[which(res$zFDM >= thresh.2)] <- 'D'
    res$df.col[which(res$zDF >= thresh.1)] <- 'C'
    res$df.col[which(res$zDF >= thresh.2)] <- 'D'
    
    res$fdm.col <- factor(res$fdm.col, levels = c('A', 'B', 'C', 'D'))
    res$df.col <- factor(res$df.col, levels = c('A', 'B', 'C', 'D'))
    
    colnames(res) <- 
      c('Scaffold', 'WindowStart', 'WindowEnd', 
        'D', 'fD', 'fDM', 'dF', 'Chromosome',
        'ID', 'windowMid', 'cumPos', 'cumStart',
        'cumEnd', 'fdm.col', 'df.col', 'zFDM',
        'zDF', 'zFDM.bonf', 'zDF.bonf', 
        'WindStart.5K', 'WindEnd.5K')
    
    res$Chromosome <- factor(res$Chromosome, 
                             levels = c(1:20,22:23))
    cols <- c('black', 'black', '#fc4e2a', '#bd0026')
    
    # Merge windows - to be plotted as vertical bands to indicate peaks of d-statistics
    if(stat == 'f_dM'){
      merged.winds <- merge.peaks(obs = res, signif.metric = 'zscore', thresh = thresh, fstat = 'fDM')
    }else{
      if(stat == 'd_f'){
        merged.winds <- merge.peaks(obs = res, signif.metric = 'zscore', thresh = thresh, fstat = 'dF')
      }else{
        merged.winds <- merge.peaks(obs = res, signif.metric = 'zscore', thresh = thresh, fstat = c('fDM', 'dF'))
      }
    }
    
    chrom.mids <- aggregate(cumPos~Chromosome, FUN=mean, data = res)
    
    chrom.begins <- c()
    chrom.ends <- c()
    for(i in 1:22){
      chr.a <- unique(as.character(res$Chromosome))[i]
      chr.b <- unique(as.character(res$Chromosome))[i+1]
      if(i == 1){
        chrom.begins[i] <- 0
        chrom.ends[i] <- min(res$cumStart[which(res$Chromosome == chr.b)])-1
      }else{
        if(i < 22){
          chrom.begins[i] <- min(res$cumStart[which(res$Chromosome == chr.a)])
          chrom.ends[i] <- min(res$cumStart[which(res$Chromosome == chr.b)])-1
        }else{
          chrom.begins[i] <- min(res$cumStart[which(res$Chromosome == chr.a)])
          chrom.ends[i] <- max(res$cumStart[which(res$Chromosome == chr.a)])
        }
      }
    }
    chroms <- 
      data.frame(chrom = c(1:20,22:23), 
                 begin = chrom.begins,
                 end = chrom.ends,
                 col = rep(c('a', 'b'), 11))
    
    chrom.cols <- c('#e0e0e0', '#d0d0d0')
    
    if(stat == 'f_dM'){
      if(length(which(res$zFDM >= 4)) > 0){
        fdm.winds <- merged.winds[which(merged.winds$stat == 'fDM'),]
        fdm <-
          ggplot() + 
          geom_rect(data = chroms, aes(xmin = begin, xmax = end, fill = col), 
                    ymin=0, ymax=Inf, alpha = 0.15) +
          scale_fill_manual(values = c('white', 'grey10')) + 
          new_scale(new_aes = 'fill') + 
          geom_vline(data = fdm.winds, aes(xintercept = CumPos),
                     alpha = 0.5, size = 1.5, color = '#DF8F443F') +
          geom_point(data = res, aes(cumPos, fDM, colour = fdm.col, 
                                     size = fdm.col, alpha = fdm.col)) +
          scale_alpha_manual(values = c(0.5, 0.5, 1, 1)) +
          scale_color_manual(values = cols) +
          scale_size_manual(values = c(0.75,0.75,1,1.25)) +
          xlab('Chromosome') + 
          labs(title = paste0(spp1, " - ", spp2)) +
          ylab(expression(F[dM])) +
          theme_classic(base_size = 24) +
          scale_x_continuous(breaks = chrom.mids$cumPos,
                             labels = chrom.mids$Chromosome, 
                             expand = c(0.025,0.025)) +
          theme(axis.text.x = element_text(size = 20), 
                axis.title.x = element_blank(), 
                axis.title.y = element_text(face = "italic"), 
                plot.title = element_text(face = "italic", size = 18),
                legend.position = 'none') 
      }else{
        fdm <-
          ggplot() + 
          geom_rect(data = chroms, aes(xmin = begin, xmax = end, fill = col), 
                    ymin=0, ymax=Inf, alpha = 0.15) +
          scale_fill_manual(values = c('white', 'grey10')) + 
          new_scale(new_aes = 'fill') + 
          geom_point(data = res, aes(cumPos, fDM, colour = fdm.col, 
                                     size = fdm.col, alpha = fdm.col)) +
          scale_alpha_manual(values = c(0.5, 0.5, 1, 1)) +
          scale_color_manual(values = cols) +
          scale_size_manual(values = c(0.75,0.75,1,1.25)) +
          xlab('Chromosome') + 
          labs(title = paste0(spp1, " - ", spp2)) +
          ylab(expression(F[dM])) +
          theme_classic(base_size = 24) +
          scale_x_continuous(breaks = chrom.mids$cumPos,
                             labels = chrom.mids$Chromosome, 
                             expand = c(0.025,0.025)) +
          theme(axis.text.x = element_text(size = 20), 
                axis.title.x = element_blank(), 
                axis.title.y = element_text(face = "italic"), 
                plot.title = element_text(face = "italic", size = 18),
                legend.position = 'none')
      }
      
      res <- list(SiteDstats = res, MergedWinds = merged.winds, 
                  fdm.plot = fdm)
    }else{
      if(length(which(res$zFDM >= 4)) > 0){
        df.winds <- merged.winds[which(merged.winds$stat == 'dF'),]
        df <-
          ggplot() + 
          geom_rect(data = chroms, aes(xmin = begin, xmax = end, fill = col), 
                    ymin=0, ymax=Inf, alpha = 0.2) +
          scale_fill_manual(values = c('white', 'grey10')) + 
          new_scale(new_aes = 'fill') + 
          geom_vline(data = df.winds, aes(xintercept = CumPos),
                     alpha = 0.5, size = 1.5, color = '#DF8F443F') +
          geom_point(data = res, aes(cumPos, dF, colour = df.col, 
                                     size = df.col, alpha = df.col)) +
          scale_alpha_manual(values = c(0.5, 0.5, 1, 1)) +
          scale_color_manual(values = cols) +
          scale_size_manual(values = c(0.75,0.75,1,1.25)) +
          xlab('Chromosome') + 
          labs(title = paste0(spp1, " - ", spp2)) +
          ylab(expression(d[F])) +
          theme_classic(base_size = 24) +
          scale_x_continuous(breaks = chrom.mids$cumPos,
                             labels = chrom.mids$Chromosome, 
                             expand = c(0.025,0.025)) +
          theme(axis.text.x = element_text(size = 20), 
                axis.title.x = element_blank(), 
                axis.title.y = element_text(face = "italic"), 
                plot.title = element_text(face = "italic", size = 18),
                legend.position = 'none')
      }else{
        df <-
          ggplot() + 
          geom_rect(data = chroms, aes(xmin = begin, xmax = end, fill = col), 
                    ymin=0, ymax=Inf, alpha = 0.2) +
          scale_fill_manual(values = c('white', 'grey10')) + 
          new_scale(new_aes = 'fill') + 
          geom_point(data = res, aes(cumPos, dF, colour = df.col, 
                                     size = df.col, alpha = df.col)) +
          scale_alpha_manual(values = c(0.5, 0.5, 1, 1)) +
          scale_color_manual(values = cols) +
          scale_size_manual(values = c(0.75,0.75,1,1.25)) +
          xlab('Chromosome') + 
          labs(title = paste0(spp1, " - ", spp2)) +
          ylab(expression(d[F])) +
          theme_classic(base_size = 24) +
          scale_x_continuous(breaks = chrom.mids$cumPos,
                             labels = chrom.mids$Chromosome, 
                             expand = c(0.025,0.025)) +
          theme(axis.text.x = element_text(size = 20), 
                axis.title.x = element_blank(), 
                axis.title.y = element_text(face = "italic"), 
                plot.title = element_text(face = "italic", size = 18),
                legend.position = 'none')
      }
      
      
      res <- list(SiteDstats = res, MergedWinds = merged.winds, 
                  df.plot = df)
    }
    return(res)
  }