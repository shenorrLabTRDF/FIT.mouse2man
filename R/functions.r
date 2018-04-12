#' Check input file format
#'
#' This function validates the format of the input mouse gene expression file is correct for process by FIT.
#' It will output an error message int he following cases:
#' (a) Less than 80% of the gene IDs are Entrez IDs  (for RNAseq data and microarray data)
#' (b) Sample names don't start with "c_" or "d_"  (for microarray daat only)
#' (c) There are at least 3 disease samples and 3 control samples  (for microarray daat only)
#' (d) The data i snot log-transformed (the range of values are either <0 or >100)  (for microarray daat only)
#' @param MouseData The mouse data
#' @param DataType 'microarray' or 'rnaseq'
CheckFormat = function(MouseData, DataType)
{
  MM_entrez = MM_Entrez_symbol_desc[,"MM.Entrez"]
  if (DataType=="microarray") names = rownames(MouseData)
  else names = MouseData$MM.Entrez
  
  per = sum(names %in% MM_entrez)*100/length(names)
  if (per<80) stop("Error: Data not in a correct format: Less than 80% of the row names are Entrez ID")
  else
  {
    if(DataType=="rnaseq") err = "\n\nSuccess: The data is in the correct format.\nThe next step is to run the analysis."
    else
    {
      samp_names = colnames(MouseData)
      if(length(grep("c_*|d_*", samp_names, perl = T, invert = T))>0)
        stop("Error: The data not in a correct format: All sample names (colnames) should start with c_ or d_.")
      else 
        if(sum(grepl("c_*", samp_names, perl = T))<3)
          stop("Error: The data not in a correct format: There are less than 3 control samples.")
        else
          if(sum(grepl("d_*", samp_names, perl = T))<3)
            stop("Error: The data not in a correct format: There are less than 3 disease samples.")
          else
            if ((range(MouseData)[1]<0) || range(MouseData)[2]>100)
              stop("Error: The data not in a correct format: It is not logged transformed.")
            else
              err = "Success: The data is in the correct format."
    }
  }
  return(err)
}



#' Pre-processing of mouse data for predictions by FIT (for microarray data only)
#' 
#' Pre-processing includes 4 steps:
#' a) Computing fold-change per gene
#' b) Computing Z-scores per gene
#' c) Computing Z-tests per gene
#' d) Merging human orthologs
#' @param MouseData The mouse data.
PreProcess = function(MouseData)
{
  dis_samp = grep("d_*", colnames(NewMouse), perl = T)
  cont_samp = grep("c_*", colnames(NewMouse), perl = T)
  FC_Mouse = apply(NewMouse, 1L, function(row) {mean(row[dis_samp], na.rm = T) - mean(row[cont_samp], na.rm = T)})
  Zscore_Mouse = apply(NewMouse, 2L, function(col) {(col - mean(col, na.rm = T))/sd(col, na.rm = T) })
    
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
  colnames(comb_data)=c("gene", "FC", "EffectSize")
  rownames(comb_data) = comb_data[,"gene"]
  comb_data = comb_data[,-1]
  comb_data = merge(comb_data, MGD_orthologs, by.x=0, by.y="Mouse", all.x=T, all.y=F)
  colnames(comb_data) = c("MM.Entrez", "FC", "EffectSize", "HS.Entrez")
  
  return(comb_data)
}



#' Compute FIT predictions
#' 
#' This function computes the predictions from a pre-processed mouse dataset and the slopes computed for the reference data. 
#' In the process confidence intervals are computed as well per gene.
#' @param NewMouse_df The pre-processed mouse dataset
#' #' @return A dataframe including the following columns:
#'   * 
ComputePredictions = function(NewMouse_df, DataType)
{
  # Computing predicitons
  predictions = sapply(NewMouse_df$MM.Entrez, function(g)
  {
    curr_slopes = slopes_per_gene_V2.0[[g]]
    curr_slopes = curr_slopes[!is.na(curr_slopes)]
    if(length(curr_slopes) != 100) curr_slopes[(length(curr_slopes)+1):100]=NA
    
    curr_MM = subset(NewMouse_df, MM.Entrez==g, "EffectSize")[,1]
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
  final = merge(final, NewMouse_df, by.x = 0, by.y="MM.Entrez")
  if(DataType=="microarray") colnames(final)[c(1,9:10)] = c("MM.Entrez", "Orig_FC", "Orig_Ztest")
  else colnames(final)[c(1,9)] = c("MM.Entrez", "Orig_Ztest")
  
  # Combining with human genes and details
  conv = HS_MM_Symbol_Entrez
  final_ann = merge(final, conv, by.x="MM.Entrez", by.y= "Mouse.Ortholog", all.x=T, all.y=F)
  if(DataType=="microarray")
  {
    colnames(final_ann) = c("Mouse.Entrez", "FIT_prediction",  "CI_low", "CI_high", "CI_size", "CI_percentile", "FIT_percentile", 
                        "UpDown", "Mouse_FoldChange", "Mouse_EffectSize", "Human.Entrez", "Human.symbol", "Human.Entrez", "Mouse.symbol")
    final_ann = final_ann[,c("Mouse.Entrez","Human.Entrez", "Mouse.symbol",  "Human.symbol",
                             "Mouse_FoldChange", "Mouse_EffectSize", "FIT_prediction","FIT_percentile",
                             "UpDown", "CI_low", "CI_high", "CI_size", "CI_percentile")]
  }
  else
  {
    colnames(final_ann) = c("Mouse.Entrez", "FIT_prediction",  "CI_low", "CI_high", "CI_size", "CI_percentile", "FIT_percentile", 
                            "UpDown", "Mouse_EffectSize", "Human.symbol", "Human.Entrez", "Mouse.symbol")
    final_ann = final_ann[,c("Mouse.Entrez","Human.Entrez", "Mouse.symbol",  "Human.symbol",
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
  if(DataType=="microarray") MouseData = read.table(MouseFile, sep=",", header=T, row.names = 1)
  else 
    {
      MouseData = read.table(MouseFile, sep=",", header=T)
      colnames(MouseData) = c("MM.Entrez","EffectSize")
      MouseData$MM.Entrez = as.character(MouseData$MM.Entrez)
    }
  Err = CheckFormat(MouseData, DataType)
  message(Err)
  
  if (DataType=="microarray") NewMouse_df = MouseData[rownames(MouseData) %in% names(slopes_per_gene_V2.0),]
  else NewMouse_df = MouseData[MouseData$MM.Entrez %in% names(slopes_per_gene_V2.0),]
  message("\nInitial number of genes: ",nrow(MouseData), "\nNumber of genes for which FIT predictions will be calculated: ", nrow(NewMouse_df))
  
  if (DataType=="microarray") 
  {
    message("\nStep 2:\nPreprocessing the input data: computing fold-change, z-scores and z-test per gene.")
    NewMouse_df = PreProcess(NewMouse_df)
  } 
    
  if (DataType=="microarray") message("\nStep 3:\nPredicting human relevant genes using the FIT model.")
  else message("\nStep 2:\nPredicting human relevant genes using the FIT model.")
  ComputePredictions(NewMouse_df, DataType)
}



#' Download Reference Data
#' 
#' This function allows access to the latest reference data used by the FIT model
#' @return A data.frame containing the following columns:
#' \itemize{
#'   \item{\strong{Disease} - Disease name.}
#'   \item{\strong{DataType} - <Technology>_<RF/ST> - \strong{Technology} can be either RNAseq or Microarrays.
#'          \strong{RF/ST} can either be \strong{RF} (Reference) denoting a cross-species pairing (CSP) where the human and mouse datasets were directly 
#'          contrasted to one another in a publication authored by the researchers who had generated at least one of the datasets, 
#'          or \strong{ST} (Standard) denoting a CSP which consisted of human datasets and corresponding datasets of a mouse model of the 
#'          human disease in a separate study.}
#'   \item{\strong{CSP_ID} - ID of the CSP.} 
#'   \item{\strong{MM.Entrez, HS.Entrez} - Gene IDs.}
#'   \item{\strong{FC.HS, FC.MM} - Mouse (MM) and human (HS) fold-change for the CSP.}
#'   \item{\strong{EffSize.HS, EffSize.MM} - Mouse and human effect sizes for the CSP. 
#'                 Z-tests for microarrays, estimated fold-change, as computed in Sleuth for RNA-seq.}
#'   \item{\strong{qval.HS, qval.MM} - Mouse (MM) and human (HS) q-values for the CSP.}
#'   }
#' @export
GetRefData = function()
{
  AllData_V2.0
}


#' Download sampel data
#' 
#' This function allows download of sample data that can be used by FIT, in CSV format.
#' The file will be saved in the name "SampleRNAseq.csv" or "SampleMicroarray.csv" depending on the DataType chosen.
#' @param DataType Either "microarray" or "rnaseq", depending on the technology by which the data was assayed.
#' @export
GetSampleData = function(DataType)
{
  if(DataType == "rnaseq") write.table(RNAseq_sample, "SampleRNAseq.csv", sep=",", quote = F, row.names = F)
  else if(DataType == "microarray") write.table(microarray_sample, "SampleMicroarray.csv", sep=",", quote = F, row.names = F)
  else  stop("Error: DataType should be 'rnaseq' or 'microarray'.")
}



#' Run FIT improvement prediction classifier 
#' 
#' The SVM classifier predicts whether FIT will be able to improve a specific mouse data
#' @param MouseFile File name that includes the mouse data, in CSV format 
#' @param qval the q-value cuttoff the user will use to interpret FIT's predictions. (default= 0.1)
#' @param FC the fold-change cuttoff the user will use to interpret FIT's predictions, given as fraction from the top. For example, 
#'            0.15 denotes the top 15\% of genes with highest fold-change. (default= 0.15)
#' @export
RunClassifier = function(MouseFile, qval=0.1, FC=0.15)
{
  # Input checks
  if(!file.exists(MouseFile)) stop(paste0("The file ",MouseFile," doesn't exist."))
  qvals = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1)
  FCs = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25,0.3,0.35,0.4)
  if(!qval %in%  qvals) stop(paste0("q-value should be one of the following:\n", paste(qvals, collapse = " ")))
  if(!FC %in%  FCs) stop(paste0("The fold-change should be one of the following:\n", paste(FCs, collapse = " ")))
  
  MouseData = read.table(MouseFile, sep=",", header=T)  # load("../data/microarray_sample.rda"); MouseData= microarray_sample; rm(microarray_sample)
  
  # Creating PC point from input mouse data
  rotations_init = pca_rotations  # load("../data/classifier/pca_rotations.rda"); rotations_init = pca_rotations; rm(pca_rotations)
  intersection_genes = rownames(MouseData)[rownames(MouseData) %in% rownames(rotations_init)]
  message("The mouse data contains ", nrow(MouseData), " genes.\nThe classifier can be based on ", nrow(rotations_init)," genes.",
          "\nThe current run will be based on ", length(intersection_genes), " genes (intersection between the current data and the classifier set of genes.")
  rotations = rotations_init[intersection_genes,]
  MouseData = MouseData[intersection_genes,]
  
  dis_samp = grep("d_*", colnames(MouseData), perl = T)
  cont_samp = grep("c_*", colnames(MouseData), perl = T)
  FC_Mouse = t(as.data.frame(apply(MouseData, 1L, function(row) {mean(row[dis_samp], na.rm = T) - mean(row[cont_samp], na.rm = T)})))
  
  MM_pca_point = FC_Mouse %*% rotations[,1:50]
  
  # Running classifier
  #require("e1071")
  best_mod = best_models  # load("../data/classifier/best_models.rda"); best_mod = best_models; rm(best_models)
  
  classifier = best_mod[[paste0(FC, "_", qval)]]
  pred_res = as.character(predict(classifier, newdata = MM_pca_point))

  if(pred_res == 0) message("It is unlikely FIT will be able to improve this dataset.")
  else message("FIT will likely improve this dataset.")
  
  message("See the performance results of the classifier to identify the performance of the classifier in the selected set of theesholds (Fold-change=",
          FC,", q-value=",qval,")")
  ShowClassifierPerformance()
  
}


#' ShowClassifierPerformance
#' 
#' Show performance of SVM classifier (as an image, taken from the paper)
#' 
#' @export
ShowClassifierPerformance = function()
{
  #require(imager)
  im <- load.image("../data/classifier/ClassifierPerformance.PNG")
  plot(im, axes = F, )
}




