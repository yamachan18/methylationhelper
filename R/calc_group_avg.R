#' Calculate the group average for a given beta matrix and a list of groupings
#'
#' @param beta A matrix with missing values, where rows are CpGs and columns are samples.
#' @param l_groupings A list where each list element is a vector of sample names corresponding to a group.
#' The list elements should contain some subset of the column names of beta.
#' 
#' @return A matrix with the average beta values for each group, where rows are CpGs and columns are groups.
#' 
#' @details This function takes a matrix with missing values and calculates the average beta values for each group.
#' For rows with all missing values, the average value will be NA.
#'
#' @examples
#' library(methylationhelper)
#' data(GSE174666)
#' l_groupings = vec_to_list(c("GSM3489491", "GSM3489492"), c("GSM3489493", "GSM3489494"))
#' 
#' @export
calc_group_avg <- function(betadf, l_groupings){
    
    l_groupavg <- list()
    for (i in 1:length(l_groupings)){
      l_groupavg[[i]] <- rowMeans(betadf[, which( colnames(betadf) %in% unlist(l_groupings[[i]]))], na.rm=TRUE)
    }
    names(l_groupavg) <- names(l_groupings)
    df_groupavg <- do.call(cbind, l_groupavg)
        
  return(df_groupavg)
  
}