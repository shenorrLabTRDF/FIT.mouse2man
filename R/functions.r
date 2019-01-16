#' Check input file 
#'
#' This function validates that the input file exists and contain a valid number of genes in the first row
#' It converts the rownames to the gene names, and for microarray data removes the first column (gene names)
#' @param MouseFile The mouse data file
#' @param DataType 'microarray' or 'rnaseq'
#' @export
CheckFile = function(MouseFile, DataType)
{
  if(!(file.exists(MouseFile))) stop(paste0("Error: input file (",MouseFile,") doesn't exist"))
  data = read.table(MouseFile, sep=",", header=1)
  
  if((ncol(data) ==2) & (DataType=="microarray")) stop("Error: It seems like you uploaded RNAseq data but picked a 'microarray' datatype.")
  if((ncol(data) >2) & (DataType=="rnaseq")) stop("Error: It seems like you uploaded microarray data but picked a 'RNAseq' datatype.")
  
  if(any(duplicated(data[,1]))) stop("Error: The mouse data contains duplicated gene names.")
  if(any(is.na(data[,1]))) stop("Error: The mouse data contains missing gene names.")
  rownames(data) = data[,1]
  if(DataType=="microarray") data = data[,-1]
  
  data
}


#' Check input file format
#'
#' This function validates the format of the input mouse gene expression file is correct for process by FIT.
#' It will output an error message int he following cases:
#' (a) Less than 80% of the gene IDs are Entrez IDs  (for RNAseq data and microarray data)
#' (b) Sample names don't start with "c_" or "d_"  (for microarray daat only)
#' (c) There are at least 3 disease samples and 3 control samples  (for microarray daat only)
#' (d) The data i snot log-transformed (the range of values are either <0 or >100)  (for microarray daat only)
#' @param MouseData The mouse data, in CSV format
#' @param DataType 'microarray' or 'rnaseq'
#' @export
CheckFormat = function(MouseData, DataType)
{
  data(MM_Entrez_symbol_desc)
  MM_entrez = MM_Entrez_symbol_desc[,"MM.Entrez"]
  MouseData_genes = rownames(MouseData)
  
    # Checking format
  per = sum(MouseData_genes %in% MM_entrez)*100/length(MouseData_genes)
  if (per<80) stop("Error: The data is not in the  correct format - less than 80% of the row names are Entrez IDs.")
  else
  {
    if(DataType=="rnaseq") return("Success: The data is in the correct format.")
    else
    {
      samp_names = colnames(MouseData)
      if(length(grep("c_*|d_*", samp_names, perl = T)) != length(samp_names))
        stop("Error: The data is not in the correct format - all sample names (colnames) should start with c_ or d_.")
      else 
        if(sum(grepl("c_*", samp_names, perl = T))<3)
          stop("Error: The data is not in the correct format - there are less than 3 control samples.")
        else
          if(sum(grepl("d_*", samp_names, perl = T))<3)
            stop("Error: The data is not in the correct format - there are less than 3 disease samples.")
          else
            if ((range(MouseData)[1]<0) || range(MouseData)[2]>100)
              stop("Error: The data is not in the correct format - the values are not logged transformed.")
    }
  }
  return("Success: The data is in the correct format.")
}



#' Pre-processing of mouse data for predictions by FIT 
#' 
#' First, the genes are filtered only to those that have slopes in the model.
#' 
#' Pre-processing for microarrays includes 4 additional steps:
#' a) Computing fold-change per gene
#' b) Computing Z-scores per gene
#' c) Computing Z-tests per gene
#' d) Merging human orthologs
#' @param MouseData The mouse data.
#' @param DataType 'microarray' or 'rnaseq'
#' @export
PreProcess = function(MouseData, DataType)
{
  data(slopes_per_gene_V2.0)
  NewMouse_df = MouseData[rownames(MouseData) %in% names(slopes_per_gene_V2.0),]
  message("\nInitial number of genes: ",nrow(MouseData), "\nNumber of genes for which FIT predictions will be calculated: ", nrow(NewMouse_df))
  
  if (DataType == "rnaseq") 
  {
    colnames(NewMouse_df)=c("MM.Entrez", "EffectSize")
    return(NewMouse_df)
  }
  else
  {
    dis_samp = grep("d_*", colnames(NewMouse_df), perl = T)
    cont_samp = grep("c_*", colnames(NewMouse_df), perl = T)
    FC_Mouse = apply(NewMouse_df, 1L, function(row) {mean(row[dis_samp], na.rm = T) - mean(row[cont_samp], na.rm = T)})
    Zscore_Mouse = apply(NewMouse_df, 2L, function(col) {(col - mean(col, na.rm = T))/sd(col, na.rm = T) })
      
    dis_n = length(dis_samp)
    con_n = length(cont_samp)
    cont_sd = apply(Zscore_Mouse[,cont_samp],1L, sd, na.rm = T)
    dis_sd = apply(Zscore_Mouse[,dis_samp],1L, sd, na.rm = T)
      
    n= nrow(Zscore_Mouse)
    Z_test = sapply(X = rownames(Zscore_Mouse),USE.NAMES = F, FUN =  function(g)
    {
      numerator = mean(Zscore_Mouse[g, dis_samp], na.rm = T) - mean(Zscore_Mouse[g, cont_samp], na.rm = T)
      denominator = sqrt(((cont_sd[g]^2)/con_n) + ((dis_sd[g]^2)/dis_n))
      numerator/denominator
    })
    
    Ztest_mean = mean(Z_test)
    Ztest_mean_sd = sd(Z_test)
    Z_test_standard = (Z_test - Ztest_mean)/Ztest_mean_sd
    
    comb_data = merge(FC_Mouse, Z_test_standard, by=0) 
    colnames(comb_data)=c("MM.Entrez", "FC", "EffectSize")
    rownames(comb_data) = comb_data$MM.Entrez
    comb_data = comb_data[,-1]
    return(comb_data)
  }
  
  
}



#' Compute FIT predictions
#' 
#' This function computes the predictions from a pre-processed mouse dataset and the slopes computed for the reference data. 
#' In the process confidence intervals are computed as well per gene.
#' @param NewMouse_df The pre-processed mouse dataset
#' @param DataType The technology used: microarray/ rnaseq
#' @return A dataframe including gene IDs, FIT's prediction, CI size, FIT percentile, mouse and human FC and effect-size 
#' @export
ComputePredictions = function(NewMouse_df, DataType)
{
  # Computing predicitons
  predictions = sapply(rownames(NewMouse_df), function(g)
  {
    curr_slopes = slopes_per_gene_V2.0[[g]]
    curr_slopes = curr_slopes[!is.na(curr_slopes)]
    if(length(curr_slopes) != 100) curr_slopes[(length(curr_slopes)+1):100]=NA
    
    curr_MM = NewMouse_df[g, "EffectSize"]
    curr_slopes + (curr_MM*curr_slopes) 
  })
  final = data.frame(mean_pred=colMeans(predictions, na.rm=T))
  
  # Computing confidence intervals
  final$Conf_low <- apply(predictions, 2, quantile, probs = 0.05, na.rm=T)
  final$Conf_high <- apply(predictions, 2, quantile, probs = 0.95, na.rm=T)
  final$CI_size = final$Conf_high - final$Conf_low
  tot=nrow(final)
  final$CI_percentile = round(rank(final$CI_size)*100/tot,2)
  
  final$pred_percentile = round(rank(abs(final$mean_pred))*100/tot,2)
  final$UpDown = sapply(final$mean_pred, function(x) ifelse(x>0,"+","-"))
  
  # Adding original data
  final = merge(final, NewMouse_df, by.x = 0, by.y=0)
  if(DataType=="microarray") colnames(final)[c(1,9:10)] = c("MM.Entrez", "Orig_FC", "Orig_Ztest")
  else 
  {
    colnames(final)[c(1,10)] = c("MM.Entrez", "Orig_Ztest")
    final = final[,-9]
  }
  
  # Combining with human genes and details
  data(HS_MM_Symbol_Entrez)
  final_ann = merge(final, HS_MM_Symbol_Entrez, by.x="MM.Entrez", by.y= "Mouse.Ortholog", all.x=T, all.y=F)
  data(MM_Entrez_symbol_desc)
  final_ann = merge(final_ann, MM_Entrez_symbol_desc, by="MM.Entrez", all.x=T, all.y=F)
  
  if(DataType=="microarray")
  {
    final_ann = final_ann[,-14]
    colnames(final_ann) = c("Mouse.Entrez", "FIT_prediction",  "CI_low", "CI_high", "CI_size", "CI_percentile", "FIT_percentile", 
                        "UpDown", "Mouse_FoldChange", "Mouse_EffectSize", "Human.symbol", "Human.Entrez", "Mouse.symbol", "Description")
    final_ann = final_ann[,c("Mouse.Entrez","Human.Entrez", "Mouse.symbol",  "Human.symbol","Description",
                             "Mouse_FoldChange", "Mouse_EffectSize", "FIT_prediction","FIT_percentile",
                             "UpDown", "CI_low", "CI_high", "CI_size", "CI_percentile")]
  }
  else
  {
    final_ann = final_ann[,-13]
    colnames(final_ann) = c("Mouse.Entrez", "FIT_prediction",  "CI_low", "CI_high", "CI_size", "CI_percentile", "FIT_percentile", 
                            "UpDown", "Mouse_EffectSize", "Human.symbol", "Human.Entrez", "Mouse.symbol", "Description")
    final_ann = final_ann[,c("Mouse.Entrez","Human.Entrez", "Mouse.symbol",  "Human.symbol","Description",
                             "Mouse_EffectSize", "FIT_prediction","FIT_percentile",
                             "UpDown", "CI_low", "CI_high", "CI_size", "CI_percentile")]
  }
  final_ann = final_ann[order(abs(final_ann$FIT_prediction), decreasing = T), ]
  final_ann
}



#' Run FIT pipeline
#' 
#' This function runs the whole FIT pipeline: checks input file format, pre-processes data (for microarray data only) 
#' and computes predictions.
#' @param MouseFile File name that includes the mouse data, in CSV format
#' @param DataType Either "microarray" or "rnaseq", depending on the technology by which the data was assayed.
#' @return A data.frame containing the following columns:
#' \itemize{
#'   \item{\strong{Mouse.Entrez, Human.Entrez, Mouse.symbol, Human.symbol} - Gene IDs}
#'   \item{\strong{Mouse_FoldChange} - Mouse fold-change as computed from the mouse input data}
#'   \item{\strong{Mouse_Ztest} - Mouse Z-test values as computed from the mouse input data}
#'   \item{\strong{FIT_prediction} - Human effect-size prediction by the FIT model}
#'   \item{\strong{FIT_percentile} - Percentiles of absolute values of FIT's predictions}
#'   \item{\strong{UpDown} - Sign of prediction}
#'   \item{\strong{CI_low, CI_low, CI_size, CI_percentile} - Confidenc einterval values (low, high), overall size (high-low) and percentile}
#'   }
#' @export
FIT = function(MouseFile, DataType)
{
  if((DataType != "microarray") & (DataType != "rnaseq")) stop("Error: DataType should be 'rnaseq' or 'microarray'.")
  if(!file.exists(MouseFile)) stop(paste0("The file ",MouseFile," doesn't exist."))
    
  message("Step 1:\nUploading input data and checking data format")
  MouseData = CheckFile(MouseFile, DataType)
  CheckFormat(MouseData, DataType)
  
  message("\nStep 2:\nPreprocessing the input data: computing fold-change, z-scores and z-test per gene.")
  NewMouse_df = PreProcess(MouseData, DataType)
    
  message("\nStep 3:\nPredicting human relevant genes using the FIT model.")
  ComputePredictions(NewMouse_df, DataType)
}





#' Run FIT improvement prediction classifier 
#' 
#' The SVM classifier predicts whether FIT will be able to improve a specific mouse data.
#' @param MouseData The pre-processed mouse dataset, or NULL in case MouseFIle ise given
#' @param MouseFile File name that includes the mouse data (log expression per gene for all disease and control sampels), in CSV format. NULL in case MouseData is given
#' @param DataType Either "microarray" or "rnaseq", depending on the technology by which the data was assayed. NULL if MouseData is given
#' @param qval the q-value cuttoff the user will use to interpret FIT's predictions. (default= 0.1)
#' @param FC the fold-change cuttoff the user will use to interpret FIT's predictions, given as fraction from the top. For example, 
#'            0.15 denotes the top 15\% of genes with highest fold-change. (default= 0.15)
#' @param verbose A logical value defining whether detialed messages should be printed.
#' @export
RunClassifier = function(MouseData=NULL, MouseFile=NULL, DataType=NULL, qval=0.1, FC=0.15, verbose=F)
{
  # Input checks
  if(!is.null(MouseFile))
    if(!file.exists(MouseFile)) stop(paste0("The file ",MouseFile," doesn't exist."))
  if((!is.null(MouseData)) & (!is.null(MouseFile)))
    stop("Error: This function should either get a MouseData data.frame, or MouseFIle character vector")
  if(is.null(DataType))
    stop("Error: Illegal value for DataType")
  if((DataType!="microarray") & (DataType!="rnaseq"))
    stop("Error: Illegal value for DataType")
  
  qvals = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1)
  FCs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25,0.3,0.35,0.4)
  if(!qval %in%  qvals) stop(paste0("q-value should be one of the following:\n", paste(qvals, collapse = " ")))
  if(!FC %in%  FCs) stop(paste0("The fold-change should be one of the following:\n", paste(FCs, collapse = " ")))
  
  if(is.null(MouseData))
  {
    MouseData = CheckFile(MouseFile, DataType)
    res = CheckFormat(MouseData , DataType)
    M1 = res 
  } else M1=""
  
  # Creating PC point from input mouse data
  data(pca_rotations)
  intersection_genes = rownames(MouseData)[rownames(MouseData) %in% rownames(pca_rotations)]
  M1 = paste0(M1, "\nThe mouse data contains ", nrow(MouseData), " genes.\nThe classifier can be based on ", nrow(pca_rotations)," genes.",
          "\nThe current run will be based on ", length(intersection_genes), " genes (intersection between the current data and the classifier set of genes.")
  rotations = pca_rotations[intersection_genes,]
  MouseData = MouseData[intersection_genes,]
  
  if(DataType == "microarray") 
  {
    dis_samp = grep("d_*", colnames(MouseData), perl = T)
    cont_samp = grep("c_*", colnames(MouseData), perl = T)
    FC_Mouse = t(as.data.frame(apply(MouseData, 1L, function(row) {mean(row[dis_samp], na.rm = T) - mean(row[cont_samp], na.rm = T)})))
  }
  else 
    FC_Mouse = as.vector(MouseData[,2])
  
  MM_pca_point = FC_Mouse %*% rotations[,1:50]
  
  # Running classifier
  data(best_models)
  classifier = best_models[[paste0(FC, "_", qval)]]
  pred_res = as.numeric(as.character(predict(classifier, newdata = MM_pca_point)))

  # Providing results
  M2 = "\n*************************  Classifier prediction  **************************\n"
  if(pred_res == 0) M2 = paste0(M2, "It is unlikely FIT will be able to improve this dataset.\n")
  else M2 = paste0(M2, "FIT will likely improve this dataset.")
  M2 = paste0(M2, "****************************************************************************\n")
  
  if(verbose)
  {
    message(M1, "\n", M2)
    message("See the performance results of the classifier to identify the performance of the classifier in the selected set of theesholds (Fold-change=",
            FC,", q-value=",qval,")")
    ShowClassifierPerformance()
  }
  else
    return(list(pred_res, M1, M2))
}


#' ShowClassifierPerformance
#' 
#' Show performance of SVM classifier (as an image, taken from the paper)
#' 
#' @export
ShowClassifierPerformance = function()
{
  im <- load.image("inst/ClassifierPerformance.PNG")
  plot(im, axes = F)
}




