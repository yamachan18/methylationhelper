#' Function for creating a list of subsetted dataframes from a dataframe of beta values based on a list of CpGs and groupings.
#' Subset the beta data frame based on the list of CpGs and groupings.
#'
#' @param betadf A data frame containing beta values for CpGs.
#' @param list_of_cpgs A list of CpG IDs to subset the beta data frame by.
#' @param list_of_groupings A list of groupings corresponding to the columns in betadf.
#' @param name_of_reference The name of the reference group to use for subsetting the betadf.
#' @return A list of subsetted data frames corresponding to each CpG and grouping combination.
#' @details This function subsets the beta data frame based on the list of CpGs and groupings, and returns a list of data frames
#' corresponding to each CpG and grouping combination. The returned list can be used for performing statistical tests per CpG.
#'
#' @export

list_of_beta_subset <- function(betadf, list_of_cpgs, list_of_groupings, name_of_reference = NULL){
  
  if (is.null(name_of_reference)){
    # query
      l_df = list()
      for (i in 1:length(list_of_cpgs)){
        l_df[[i]] = betadf[rownames(betadf) %in% names(list_of_cpgs[[i]]), 
                                colnames(betadf) %in% list_of_groupings[[i]]]
      }
      names(l_df) <- names(list_of_groupings)
      return(l_df)
  }else{
    # reference
      w = which(names(list_of_groupings) %in% name_of_reference)
      
      l_df = list()
      for (i in 1:length(list_of_cpgs)){
        l_df[[i]] = betadf[rownames(betadf) %in% names(list_of_cpgs[[i]]), 
                                colnames(betadf) %in% list_of_groupings[[w]]] # note this is different than above
      }
      names(l_df) <- rep(names(list_of_groupings)[w], length(l_df)) # different
      return(l_df)
    
  }

}