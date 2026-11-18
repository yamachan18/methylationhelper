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


# this plots for all colnames in query_minus_immature
scatter_based_probelist <- function(query_minus_immature, reference_averages, cutoff_value=0.2, whichslice){
  
  cutoff <- cutoff_value # can change this part manually
  margin <- cutoff * sqrt(2)
  
  if (whichslice == "hypo"){
    # if colorslice = TRUE
    # default is "hypomethylation" or bottom center slice
    # default color is blue
    # select big bottom diagonal first, then difference of the differences: (diff TumorAvg - NBC) - (diff Mature - Immature)
    
    diffdiff<- apply(query_minus_immature, 2, function(x) x - reference_averages$mature_minus_immature)
    
    listofdfs <- list()
    
    for (i in 1:ncol(query_minus_immature)){
      # bottom diagonal: 
      filter1<-rownames(diffdiff) [which(diffdiff[,i] < -margin)]
      
      # then take the bottom horizontal
      filter2 <- rownames(query_minus_immature) [which(query_minus_immature[,i] < -cutoff)]
      
      # then subset on the vertical lines
      filter3 <- rownames(reference_averages)[which(reference_averages$mature_minus_immature > -cutoff &
                                             reference_averages$mature_minus_immature < cutoff)]
      
      df1 <- as.data.frame( intersect(filter1, intersect(filter2, filter3)) )
      colnames(df1) <- "TargetID"
      listofdfs[[i]]<-df1
      
      # rename my list with samplenames, except take out all the spaces and special characters s.t. it doesn't mess with filenaming and bash script
      names(listofdfs)[i] <- gsub("[^[:alnum:]]", "", colnames(query_minus_immature)[i])
    }
    
    return(listofdfs)

    } else if (whichslice == "hyper"){
    
      diffdiff<- apply(query_minus_immature, 2, function(x) x - reference_averages$mature_minus_immature)
    
      listofdfs <- list()
      
      for (i in 1:ncol(query_minus_immature)){
        # bottom diagonal: 
        filter1<-rownames(diffdiff) [which(diffdiff[,i] > margin)]
        
        # then take the bottom horizontal
        filter2 <- rownames(query_minus_immature) [which(query_minus_immature[,i] > cutoff)]
        
        # then subset on the vertical lines
        filter3 <- rownames(reference_averages)[which(reference_averages$mature_minus_immature > -cutoff &
                                               reference_averages$mature_minus_immature < cutoff)]
        
        df1 <- as.data.frame( intersect(filter1, intersect(filter2, filter3)) )
        colnames(df1) <- "TargetID"
        listofdfs[[i]]<-df1
        
        # rename my list with samplenames, except take out all the spaces and special characters s.t. it doesn't mess with filenaming and bash script
        names(listofdfs)[i] <- gsub("[^[:alnum:]]", "", colnames(query_minus_immature)[i])
      
      } 
      
      return(listofdfs)
      
      }else{
      #nothing for now
    }
  
  
}
