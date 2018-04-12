#' All the training data FIT is based on
#' 
#' Comprised of gene level expression data of 170 cross-species pairs (CSP).
#' 
#' @format A data frame with 1,737,091 rows and 11 columns:
#' \describe{
#'   \item{Disease}{The disease of the dataset}
#'   \item{DataType}{Details the technology(RNAseq/ microarrays) and the type of the CSP (standard dataset (ST) or reference (RF)), in the format <technology>_<type>}
#'   \item{CSP_ID}{A unique identifier for the CSP}
#'   \item{MM.Entrez}{Mouse Entrez ID}
#'   \item{HS.Entrez}{Human Entrez ID}
#'   \item{FC.HS}{Human fold-change}
#'   \item{MM.HS}{Mouse fold-change}
#'   \item{EffSize.HS}{Human effect size}
#'   \item{EffSize.MM}{Mouse effect size}
#'   \item{qval.HS}{Human q-value}
#'   \item{qval.MM}{Mouse q-value}
#' }
"AllData_V2.0"