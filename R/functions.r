#' Check input file format
#'
#' This function validates the format of the input mouse gene expression file is correct for process by FIT.
#' It will output an error message in 3 cases: 
#' (a) Less than 80% of the gene IDs are Entrez IDs; 
#' (b) Sample names don't start with "c_" or "d_"
#' (c) There are at least 3 disease samples and 3 control samples
#' (d) The data i snot log-transformed (the range of values are either <0 or >100)
#' @param MouseData The mouse data
CheckFormat = function(MouseData)
{
  conv = MM_Entrez_symbol_desc
  MM_entrez = conv[,"MM.Entrez"]
  names = rownames(MouseData)
  
  per = sum(names %in% MM_entrez)*100/length(names)
  if (per<80) return("Error: Data not in a correct format: Less than 80% of the row names are Entrez ID")
  else
  {
    samp_names = colnames(MouseData)
    if(length(grep("c_*|d_*", samp_names, perl = T, invert = T))>0)
      err="Error: The data not in a correct format: All sample names (colnames) should start with c_ or d_."
    else 
      if(sum(grepl("c_*", samp_names, perl = T))<3)
        err="Error: The data not in a correct format: There are less than 3 control samples."
      else
        if(sum(grepl("d_*", samp_names, perl = T))<3)
          err = "Error: The data not in a correct format: There are less than 3 disease samples."
        else
          if ((range(MouseData)[1]<0) || range(MouseData)[2]>100)
            err="Error: The data not in a correct format: It is not logged transformed."
          else
            err = "Success: The data is in the correct format."
  }
  return(err)
}



#' Pre-processing of mouse data for predictions by FIT
#' 
#' Pre-processing includes 4 steps:
#' a) Computing fold-change per gene
#' b) Computing Z-scores per gene
#' c) Computing Z-tests per gene
#' d) Merging human orthologs
#' @param MouseData The mouse data.
PreProcess = function(MouseData)
{
  slopes_per_gene = slopes_per_gene_V2.0
  NewMouse = MouseData[rownames(MouseData) %in% names(slopes_per_gene),]
    
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
  colnames(comb_data)=c("gene", "FC", "Ztest")
  rownames(comb_data) = comb_data[,"gene"]
  comb_data = comb_data[,-1]
  conv = MGD_orthologs
  comb_data = merge(comb_data, conv, by.x=0, by.y="Mouse", all.x=T, all.y=F)
  colnames(comb_data) = c("MM.Entrez", "FC", "Ztest", "HS.Entrez")
  
  return(comb_data)
}



#' Compute FIT predictions
#' 
#' This function computes the predictions from a pre-processed mouse dataset and the slopes computed for the reference data. 
#' In the process confidence intervals are computed as well per gene.
#' @param NewMouse_df The pre-processed mouse dataset
#' #' @return A dataframe including the following columns:
#'   * 
ComputePredictions = function(NewMouse_df)
{
  slopes = slopes_per_gene_V2.0
  
  # Computing predicitons
  predictions = sapply(NewMouse_df$MM.Entrez, function(g)
  {
    curr_slopes = slopes[[g]]
    curr_slopes = curr_slopes[!is.na(curr_slopes)]
    if(length(curr_slopes) != 100) curr_slopes[(length(curr_slopes)+1):100]=NA
    
    curr_MM = subset(NewMouse_df, MM.Entrez==g, Ztest)[,1]
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
  colnames(final)[c(1,9:10)] = c("MM.Entrez", "Orig_FC", "Orig_Ztest")
  
  # Combining with human genes and details
  conv = HS_MM_Symbol_Entrez
  final_ann = merge(final, conv, by.x="MM.Entrez", by.y= "Mouse.Ortholog", all.x=T, all.y=F)
  colnames(final_ann) = c("Mouse.Entrez", "FIT_prediction",  "CI_low", "CI_high", "CI_size", "CI_percentile", "FIT_percentile", 
                      "UpDown", "Mouse_FoldChange", "Mouse_Ztest", "Human.Entrez", "Human.symbol", "Human.Entrez", "Mouse.symbol")
  final_ann = final_ann[,c("Mouse.Entrez","Human.Entrez", "Mouse.symbol",  "Human.symbol",
                           "Mouse_FoldChange", "Mouse_Ztest", "FIT_prediction","FIT_percentile",
                           "UpDown", "CI_low", "CI_high", "CI_size", "CI_percentile")]
  final_ann = final_ann[order(abs(final_ann$FIT_prediction), decreasing = T), ]
  final_ann
}



#' Run FIT pipeline
#' 
#' This function runs the whole FIT pipeline: checks input file format, pre-processes data and computes predictions.
#' @param MouseFile File name that includes the mouse data, in CSV format
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
FIT = function(MouseFile)
{
  message("Step 1:\nUploading input data and checking data format")
  MouseData = read.table(MouseFile, sep=",", header=T, row.names = 1)
  Err = CheckFormat(MouseData)
  message(Err)
  
  message("\nStep 2:\nPreprocessing the input data: computing fold-change, z-scores and z-test per gene.")
  NewMouse_df = PreProcess(MouseData)
    
  message("\nStep 3:\nPredicting human relevant genes using the FIT model.")
  ComputePredictions(NewMouse_df)
}


