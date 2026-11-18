#' Plot scatterplots for methylation differences between query and reference samples
#'
#' @param query_minus_immature A matrix of beta values, where rows are CpGs and columns are samples.
#' @param reference_averages A data frame with the average methylation values for each CpG in the reference groups.
#' @param cutoff_value Values between 0 and 1 for plotting the cutoff lines. CpGs within the cutoff will be colored differently. Default is 0.2.
#' @param alpha_value Values between 0 and 1 controlling the transparency for the scatterplot. Default is 0.1. 1 is fully opaque.
#' @param name_immature The name of the immature group to be used on the plot axes.
#' @param name_mature The name of the mature group to be used on the plot axes.
#' @param name_tumor The name of the tumor group to be used on the plot axes.
#' @param colorslice Boolean value whether or not to color the points within the slice of interest. Defaulted to FALSE.
#' @param whichslice The slice to color ("hypo" to color the bottom center slice of aberrantly hypomethylated CpGs or 
#' "hyper" to color the top center slice of aberrantly hypermethylated CpGs).
#' @return A list of ggplot objects.
#' 
#' @details This function generates scatterplots for each sample (column) in the query_minus_immature matrix, 
#' comparing the methylation differences between the query and reference samples. The function allows for 
#' customization of cutoff values, transparency, and labeling of the axes. Additionally, it can color points within 
#' a specified slice of interest (hypomethylated or hypermethylated CpGs) based on the provided parameters.
#' 
#' @export
#' 
methylation_scatter<- function(query_minus_immature, reference_averages, cutoff_value=0.2, alpha_value = 0.1,
                               name_immature = "Immature", name_mature = "Mature", name_tumor = "Tumor", 
                               colorslice = FALSE, whichslice=NULL){
  
  cutoff <- cutoff_value # can change this part manually
  margin <- cutoff * sqrt(2)
  
  
  if (colorslice == FALSE){
  for (i in 1:ncol(query_minus_immature)){
    
    mydf <- data.frame(mature_minus_immature = reference_averages$mature_minus_immature, 
                       query_minus_immature = query_minus_immature[,i])
    mydf_narm <- tidyr::drop_na(mydf)
  
    gg<-ggplot2::ggplot(mydf_narm, aes(x=mature_minus_immature, y=query_minus_immature)) + 
    geom_scattermore(alpha= alpha_value)+ #alpha changes the density; scattemore speeds up plotting
    geom_abline(intercept=0, slope=1, linetype="dashed",color="red")+ #y=x line
    geom_abline(intercept=margin, slope=1, linetype="solid",color="red")+
    geom_abline(intercept=-margin, slope=1, linetype="solid",color="red")+
    geom_vline(xintercept = cutoff, color="red")+
    geom_vline(xintercept = -cutoff, color="red")+
    geom_hline(yintercept = cutoff, color="red")+
    geom_hline(yintercept = -cutoff, color="red")+
    xlim(-1,1)+
    ylim(-1,1)+
    theme_bw()+
    theme(aspect.ratio = 1)+ # make the plot square
    ggtitle(paste0(colnames(query_minus_immature)[i]))+
    labs(
      x= paste0(name_immature ," -> ", name_mature),
      y=paste0(name_immature ," -> ", name_tumor)
    )
  
    print(gg)
    }
  } else if (colorslice == TRUE & whichslice == "hypo"){
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
    
    
    for (i in 1:ncol(query_minus_immature)){
      colorthese<-listofdfs[[i]]$TargetID
      yeet<-ifelse(rownames(query_minus_immature) %in% colorthese, "subset", "ignore")
      mydf <- data.frame(mature_minus_immature = reference_averages$mature_minus_immature, 
                       query_minus_immature = query_minus_immature[,i], 
                       colors = yeet)
      mydf_narm <- tidyr::drop_na(mydf)
      
      
      gg<-ggplot2::ggplot(mydf_narm, aes(x=mature_minus_immature, y=query_minus_immature, col=colors)) + 
      geom_scattermore(alpha= alpha_value)+ #alpha changes the density; scattemore speeds up plotting
      geom_abline(intercept=0, slope=1, linetype="dashed",color="red")+ #y=x line
      geom_abline(intercept=margin, slope=1, linetype="solid",color="red")+
      geom_abline(intercept=-margin, slope=1, linetype="solid",color="red")+
      geom_vline(xintercept = cutoff, color="red")+
      geom_vline(xintercept = -cutoff, color="red")+
      geom_hline(yintercept = cutoff, color="red")+
      geom_hline(yintercept = -cutoff, color="red")+
      xlim(-1,1)+
      ylim(-1,1)+
      theme_bw()+
      theme(aspect.ratio = 1)+ # make the plot square
      ggtitle(paste0(colnames(query_minus_immature)[i]))+
      labs(
        x= paste0(name_immature ," -> ", name_mature),
        y=paste0(name_immature ," -> ", name_tumor)
      )
    
      print(gg)
      
    }
    
    }else if (colorslice == TRUE & whichslice == "hyper"){
    
      
      
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
    
    
    for (i in 1:ncol(query_minus_immature)){
      colorthese<-listofdfs[[i]]$TargetID
      yeet<-ifelse(rownames(query_minus_immature) %in% colorthese, "subset", "ignore")
      mydf <- data.frame(mature_minus_immature = reference_averages$mature_minus_immature, 
                       query_minus_immature = query_minus_immature[,i], 
                       colors = yeet)
      mydf_narm <- tidyr::drop_na(mydf)
      
      
      gg<-ggplot2::ggplot(mydf_narm, aes(x=mature_minus_immature, y=query_minus_immature, col=colors)) + 
      geom_scattermore(alpha= alpha_value)+ #alpha changes the density; scattemore speeds up plotting
      geom_abline(intercept=0, slope=1, linetype="dashed",color="red")+ #y=x line
      geom_abline(intercept=margin, slope=1, linetype="solid",color="red")+
      geom_abline(intercept=-margin, slope=1, linetype="solid",color="red")+
      geom_vline(xintercept = cutoff, color="red")+
      geom_vline(xintercept = -cutoff, color="red")+
      geom_hline(yintercept = cutoff, color="red")+
      geom_hline(yintercept = -cutoff, color="red")+
      xlim(-1,1)+
      ylim(-1,1)+
      theme_bw()+
      theme(aspect.ratio = 1)+ # make the plot square
      ggtitle(paste0(colnames(query_minus_immature)[i]))+
      labs(
        x= paste0(name_immature ," -> ", name_mature),
        y=paste0(name_immature ," -> ", name_tumor)
      )
    
      print(gg)
      
      }
      
      
    }else{
      #nothing for now
    }
  
  
}