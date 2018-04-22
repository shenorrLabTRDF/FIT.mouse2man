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
#'   \item{FC.MM}{Mouse fold-change}
#'   \item{EffSize.HS}{Human effect size}
#'   \item{EffSize.MM}{Mouse effect size}
#'   \item{qval.HS}{Human q-value}
#'   \item{qval.MM}{Mouse q-value}
#' }
"AllData_V2.0"

#' A data.frame containing the Entrez IDs of human and mouse orthologs
"HS_MM_Symbol_Entrez"

#' A data.frame containing the Entrez IDs of human and mouse orthologs
"MGD_orthologs"

#' A data.frame containing the symbols of human and mouse orthologs and a description per gene
"MM_Entrez_symbol_desc"

#' A list of the best predictive model for FIT's perforemance per threhsolds pair (fold-change and q-value)
"best_models"

#' The rotations needed to convert a vector of fold-changes to a dot on the PC space of all cross-species pairs in the training data
"pca_rotations"

#' The slopes derived from the linear regression FIT computed from human and mouse data, per gene. This is used for the prediciton of human effect size values.
"slopes_per_gene_V2.0"

#' Microarray sample data
"microarray_sample"

#' RNAseq sample data
"RNAseq_sample"

