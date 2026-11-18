mean_imputation = function(beta){
  
  k <- which(is.na(beta), arr.ind=TRUE)
  beta[k] <- rowMeans(beta, na.rm=TRUE)[k[,1]]
  
  return(beta)
}