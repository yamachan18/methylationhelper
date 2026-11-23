#' Convert a data frame to a list of named vectors, one for each column.
#'
#' @param df A data frame to be converted to a list of named vectors.
#' @return A list of named vectors where the name of each vector corresponds to the column name in the original data frame.
#' @details This function converts a data frame to a list of named vectors, where each vector corresponds to a column in the original data frame. The names of the vectors are taken from the row names of the data frame.
#'
#' @export

df_to_list <- function(df){
  
  l_u <- list()
  for (i in 1:ncol(df)){
    u = colnames(df)[i]
    
    # named vector
    named_vec = df[,i]
    names(named_vec) <- rownames(df)
    
    l_u[[i]] <- named_vec
    
  }
  names(l_u) <- colnames(df)
  return(l_u)
}