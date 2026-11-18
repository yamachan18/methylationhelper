#' Generate a list of CpGs based on the slice used in methylation_scatter
#'
#' @param query_minus_immature A dataframe of methylation values for the query samples minus the reference samples in the immature group.
#' @param reference_averages A data frame with the average methylation values for each CpG in the reference groups.
#' @param cutoff_value A numeric value between 0 and 1 for plotting the cutoff lines. CpGs within the cutoff will be colored differently. Default is 0.2.
#' @param whichslice The slice to color ("hypo" to color the bottom center slice of aberrantly hypomethylated CpGs or 
#' "hyper" to color the top center slice of aberrantly hypermethylated CpGs).
#' 
#' @return A list where each list element contains a dataframe of CpGs (TargetIDs) that fall within the specified slice for each sample (column) in query_minus_immature.
#' 
#' @details This function identifies CpGs that fall within a specified slice of interest (hypomethylated or hypermethylated) 
#' based on the differences between query and reference methylation values. It returns a list of dataframes, each containing 
#' the TargetIDs of CpGs that meet the criteria for each sample in the query_minus_immature matrix. The function allows for 
#' customization of cutoff values and selection of the slice of interest.
#'
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
