vec_to_list <- function(vec_groups, vec_values){
  
  l_u <- list()
  for (i in 1:length(unique(vec_groups))){
    u = unique(vec_groups)[i]
    
    l_u[[i]] <- vec_values[vec_groups %in% u]
    
  }
  names(l_u) <- unique(vec_groups)
  return(l_u)
}