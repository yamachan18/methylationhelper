#' Turn two vectors into a list
#'
#' @param vec_groups A vector of groupings, where each element corresponds to a value in vec_values.
#' @param vec_values A vector of values, where each element corresponds to a group in vec_groups.
#' @return A list with elements corresponding to each unique group in vec_groups.
#' 
#' @details This function takes two vectors, one for groupings and one for values, 
#' and returns a list where each list element corresponds to a group and each group contains values
#' from vec_values belonging to that group.
#'
#' @examples
#' mylist = vec_to_list(vec_groups=c("A","A","B","B"), vec_values=c(1,2,3,4))
#' 
#' @export


vec_to_list <- function(vec_groups, vec_values){
  
  l_u <- list()
  for (i in 1:length(unique(vec_groups))){
    u = unique(vec_groups)[i]
    
    l_u[[i]] <- vec_values[vec_groups %in% u]
    
  }
  names(l_u) <- unique(vec_groups)
  return(l_u)
}