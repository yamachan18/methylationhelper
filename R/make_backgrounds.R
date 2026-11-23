#' Generate a restricted universe of backgrounds to be used in HOMER motif enrichment analysis
#'
#' @param listofdfs A list of data frames with the CpG IDs to be used as the foreground in motif enrichment analysis.
#' @param reference A data frame with the CpG IDs for the reference set.
#' @param seed A numeric value for setting the random seed, for reproducibility.
#' @param write A logical value indicating whether to write the output to files.
#' @param write_dir The directory to write the output files if write is TRUE.
#' 
#' @return If write is TRUE, the function writes the output to files in the specified write_dir.
#' Otherwise, the function returns the generated foreground and background data frames in BED format as a list.
#' @details This function generates a set of background CpGs for each set of foreground CpGs in the list, 
#' ensuring that the background CpGs are chosen from the set difference between the full 'universe' of CpGs and the foreground CpGs.
#' We bin these background CpGs according to their rounded beta value between 0.1 and 0.9, 
#' and choose an evenly represented subset of them to represent the background.
#' The output can be written to files or returned as a list of data frames in BED format.
#' 
#'
#' @export 


make_backgrounds <- function(listofdfs, reference, seed, write = TRUE, write_dir=NULL){

l_l_HOMER <- list()  
for (i in 1:length(listofdfs)){

  Foreground <- listofdfs[[i]]
  
  ForegroundMean <- merge(reference, Foreground, by = "TargetID")
  
  BackgroundProbes <- data.frame(setdiff(reference$TargetID, Foreground$TargetID))
  colnames(BackgroundProbes) <- "TargetID"
  
  BackgroundMean <- merge(reference, BackgroundProbes, by ="TargetID")
  
  ForegroundMeanSorted <- ForegroundMean[order(ForegroundMean$refMean, decreasing = TRUE),]
  ForegroundMeanSorted$Bracket = round(ForegroundMeanSorted$refMean,1)
  NumberForeground = nrow(ForegroundMeanSorted)
  print(NumberForeground)
  ForegroundBracket10 = sum(with(ForegroundMeanSorted, Bracket == 1))/NumberForeground
  ForegroundBracket9 = sum(with(ForegroundMeanSorted, Bracket == .9))/NumberForeground
  ForegroundBracket8 = sum(with(ForegroundMeanSorted, Bracket == .8))/NumberForeground
  ForegroundBracket7 = sum(with(ForegroundMeanSorted, Bracket == .7))/NumberForeground
  ForegroundBracket6 = sum(with(ForegroundMeanSorted, Bracket == .6))/NumberForeground
  ForegroundBracket5 = sum(with(ForegroundMeanSorted, Bracket == .5))/NumberForeground
  ForegroundBracket4 = sum(with(ForegroundMeanSorted, Bracket == .4))/NumberForeground
  ForegroundBracket3 = sum(with(ForegroundMeanSorted, Bracket == .3))/NumberForeground
  ForegroundBracket2 = sum(with(ForegroundMeanSorted, Bracket == .2))/NumberForeground
  ForegroundBracket1 = sum(with(ForegroundMeanSorted, Bracket == .1))/NumberForeground
  ForegroundBracket0 = sum(with(ForegroundMeanSorted, Bracket == 0))/NumberForeground
  
  BackgroundMeanSorted <- BackgroundMean[order(BackgroundMean$refMean, decreasing = TRUE),]
  BackgroundMeanSorted$Bracket = round(BackgroundMeanSorted$refMean,1) 
  
  set.seed(seed)
  BackgroundMeanSorted$RND <- sample(1000, size = nrow(BackgroundMeanSorted), replace = TRUE) # change to make sure it's reproducible
  
  BackgroundBracket10 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == 1,])
  BackgroundBracket9 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == .9,])
  BackgroundBracket8 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == .8,])
  BackgroundBracket7 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == .7,])
  BackgroundBracket6 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == .6,])
  BackgroundBracket5 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == .5,])
  BackgroundBracket4 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == .4,])
  BackgroundBracket3 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == .3,])
  BackgroundBracket2 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == .2,])
  BackgroundBracket1 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == .1,])
  BackgroundBracket0 <- data.frame(BackgroundMeanSorted[BackgroundMeanSorted$Bracket == 0,])
  
  BackgroundBracket10 <- BackgroundBracket10[order(BackgroundBracket10$RND, decreasing = TRUE),]
  BackgroundBracket9 <- BackgroundBracket9[order(BackgroundBracket9$RND, decreasing = TRUE),]
  BackgroundBracket8 <- BackgroundBracket8[order(BackgroundBracket8$RND, decreasing = TRUE),]
  BackgroundBracket7 <- BackgroundBracket7[order(BackgroundBracket7$RND, decreasing = TRUE),]
  BackgroundBracket6 <- BackgroundBracket6[order(BackgroundBracket6$RND, decreasing = TRUE),]
  BackgroundBracket5 <- BackgroundBracket5[order(BackgroundBracket5$RND, decreasing = TRUE),]
  BackgroundBracket4 <- BackgroundBracket4[order(BackgroundBracket4$RND, decreasing = TRUE),]
  BackgroundBracket3 <- BackgroundBracket3[order(BackgroundBracket3$RND, decreasing = TRUE),]
  BackgroundBracket2 <- BackgroundBracket2[order(BackgroundBracket2$RND, decreasing = TRUE),]
  BackgroundBracket1 <- BackgroundBracket1[order(BackgroundBracket1$RND, decreasing = TRUE),]
  BackgroundBracket0 <- BackgroundBracket0[order(BackgroundBracket0$RND, decreasing = TRUE),]
  
  BackgroundNumber <- 100000
  
  while (ForegroundBracket10*BackgroundNumber > nrow(BackgroundBracket10))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  while (ForegroundBracket9*BackgroundNumber > nrow(BackgroundBracket9))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  while (ForegroundBracket8*BackgroundNumber > nrow(BackgroundBracket8))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  while (ForegroundBracket7*BackgroundNumber > nrow(BackgroundBracket7))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  while (ForegroundBracket6*BackgroundNumber > nrow(BackgroundBracket6))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  while (ForegroundBracket5*BackgroundNumber > nrow(BackgroundBracket5))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  while (ForegroundBracket4*BackgroundNumber > nrow(BackgroundBracket4))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  while (ForegroundBracket3*BackgroundNumber > nrow(BackgroundBracket3))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  while (ForegroundBracket2*BackgroundNumber > nrow(BackgroundBracket2))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  while (ForegroundBracket1*BackgroundNumber > nrow(BackgroundBracket1))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  while (ForegroundBracket0*BackgroundNumber > nrow(BackgroundBracket0))
  {
    BackgroundNumber = BackgroundNumber - 10000
  }
  
  print(BackgroundNumber)
  
  BackgroundProbesReduced <- rbind(BackgroundBracket10[1:trunc(ForegroundBracket10*BackgroundNumber),], 
                                   BackgroundBracket9[1:trunc(ForegroundBracket9*BackgroundNumber),],
                                   BackgroundBracket8[1:trunc(ForegroundBracket8*BackgroundNumber),],
                                   BackgroundBracket7[1:trunc(ForegroundBracket7*BackgroundNumber),],
                                   BackgroundBracket6[1:trunc(ForegroundBracket6*BackgroundNumber),],
                                   BackgroundBracket5[1:trunc(ForegroundBracket5*BackgroundNumber),],
                                   BackgroundBracket4[1:trunc(ForegroundBracket4*BackgroundNumber),],
                                   BackgroundBracket3[1:trunc(ForegroundBracket3*BackgroundNumber),],
                                   BackgroundBracket2[1:trunc(ForegroundBracket2*BackgroundNumber),],
                                   BackgroundBracket1[1:trunc(ForegroundBracket1*BackgroundNumber),],
                                   BackgroundBracket0[1:trunc(ForegroundBracket0*BackgroundNumber),])
  
  BackgroundProbesReduced <- BackgroundProbesReduced[order(BackgroundProbesReduced$refMean, decreasing = TRUE),]
  
  BackgroundProbesHomer <- BackgroundProbesReduced[,c("CHR","MAPINFO","END")]
  rownames(BackgroundProbesHomer) <- BackgroundProbesReduced$TargetID

  ForegroundProbesHomer <- ForegroundMeanSorted[,c("CHR","MAPINFO","END")]
  rownames(ForegroundProbesHomer) <- ForegroundMeanSorted$TargetID

  plot(ForegroundMeanSorted$refMean)
  plot(BackgroundProbesReduced$refMean)
  
  l_HOMER = list(
    "BackgroundProbesHomer" = BackgroundProbesHomer,
    "ForegroundProbesHomer" = ForegroundProbesHomer
  )

  
  if (write == TRUE & !is.null(write_dir)){
    
      subfolder_name <- names(listofdfs)[i]

        if (!dir.exists(write_dir)){
            dir.create(write_dir, recursive = TRUE)
        }

      dir.create(path=file.path(write_dir, subfolder_name), recursive=TRUE)
      
      write.csv(BackgroundProbesReduced, file = file.path(write_dir, subfolder_name, "BackgroundProbes.csv"))
      write.csv(ForegroundMeanSorted, file = file.path(write_dir, subfolder_name, "ForegroundProbes.csv"))
      write.table(BackgroundProbesHomer, file=file.path(write_dir, subfolder_name, "bkg.txt"),
                  sep="\t", col.names = F, row.names = F, quote = FALSE)
      write.table(ForegroundProbesHomer, file=file.path(write_dir, subfolder_name, "target.txt"),
                  sep="\t", col.names = F, row.names = F, quote = FALSE)
    
  }else if(write == TRUE & is.null(write_dir)){
    
    stop("Please indicate what directory to write to.")
    
  }else{
    
    l_l_HOMER[[i]] <- l_HOMER

  }
  
}

if (length(l_l_HOMER) > 0){

  names(l_l_HOMER) <- names(listofdfs)

  return(l_l_HOMER)
}

}
