#' Perform Wilcoxon rank-sum test for each CpG between reference and query groups.
#'
#' @param list_of_cpgs A list of CpG IDs for which to perform the Wilcoxon rank-sum test.
#' @param list_of_dfs_ref A list of data frames corresponding to the reference groups for each CpG.
#' @param list_of_dfs_query A list of data frames corresponding to the query groups for each CpG.
#' @param l_groupings A list where each list element is a vector of sample names corresponding to a group.
#' 
#' @return A list of data frames containing the Wilcoxon rank-sum p-values (raw, FDR-adjusted, and Bonferroni-corrected) for each CpG.
#' 
#' @details This function performs the Wilcoxon rank-sum test for each CpG between the reference and query groups.
#' The results include raw p-values, FDR-adjusted p-values, and Bonferroni-adjusted p-values for each CpG.
#'
#' @export

fun_wilcoxon = function(list_of_cpgs, list_of_dfs_ref, list_of_dfs_query){
  
  # Maybe insert some message here "evaluating XXXX samples... might take a while"
  
  l_df <- list()
  for (j in 1:length(list_of_cpgs)){
  
  cpgids = list_of_cpgs[[j]]
  mat_ref = list_of_dfs_ref[[j]]
  mat_query = list_of_dfs_query[[j]]
  
  l_wil <- list()
    for (i in 1:length(cpgids)){
      wil = suppressWarnings(wilcox.test(x= as.numeric(mat_query[i,]) , y= as.numeric(mat_ref[i,]), alternative = "two.sided")$p.value)
      l_wil[[i]]<-wil
    }
  
    pvaldf = data.frame(
     cpgid = cpgids,
     pval = unlist(l_wil),
     pval_fdr = p.adjust(unlist(l_wil), method = "fdr"),
     pval_bonferroni = p.adjust(unlist(l_wil), method = "bonferroni")
   )
    
    l_df[[j]] <- pvaldf
  
  }
  
  return(l_df)
}