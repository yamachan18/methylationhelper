#' Impute with the row mean for a matrix with missing values
#'
#' @param beta A matrix with missing values, where rows are CpGs and columns are samples.
#' @return A matrix with missing values for each row (CpG) imputed with the mean of the samples.
#' 
#' @details This function takes a matrix with missing values and imputes the missing values for each row (CpG) with the mean of the samples.
#' For rows with all missing values, the imputed value will be NA.
#'
#' @examples
#' m = matrix(1:12, nrow=4, ncol=3)
#' m[2,2] <- NA
#' m
#' m_imputed <- mean_imputation(beta=m)
#' m_imputed
#' 
#' @export

mean_imputation = function(beta){
  
  k <- which(is.na(beta), arr.ind=TRUE)
  beta[k] <- rowMeans(beta, na.rm=TRUE)[k[,1]]
  
  return(beta)
}