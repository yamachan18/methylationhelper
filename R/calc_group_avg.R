calc_group_avg <- function(betadf, l_groupings){
    
    l_groupavg <- list()
    for (i in 1:length(l_groupings)){
      l_groupavg[[i]] <- rowMeans(betadf[, which( colnames(betadf) %in% unlist(l_groupings[[i]]))], na.rm=TRUE)
    }
    names(l_groupavg) <- names(l_groupings)
    df_groupavg <- do.call(cbind, l_groupavg)
        
  return(df_groupavg)
  
}