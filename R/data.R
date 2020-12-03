#' @title Bladder Cancer Clinical Data
#' @description the clinical information (including: time, status...) of bladder cancer
#' @name blca_clinical
#' @usage blca_clinical
#' @format
#' \describe{
#'      \item{ids}{sample ID}
#'      \item{status}{sample survival status (1 = dead,0 = alive)}
#'      \item{time}{survival time (days)}
#'      \item{Gender}{the sample gender}
#'      \item{...}{other clinical informations}
#' }
NULL

#' @title Bladder Cancer gene expression profile
#' @description the gene expression value matrix of bladder cancer
#' @name blca_exp
#' @usage blca_exp
#' @format row are genes, column are samples
#'
NULL

#' @title Hypoxia-related gene set
#' @description a gene set contains the hypoxia-related genes
#' @name hypoxia_genes
#' @usage hypoxia_genes
#' @format a character vector
#'
NULL

#' @title A special case of data
#' @description a special case of data run \link{risk_score}() function
#' @name meta_dat
#' @usage meta_dat
#' @format
#' \describe{
#'      \item{ids}{sample ID}
#'      \item{time}{survival time (days)}
#'      \item{status}{sample survival status (1 = dead,0 = alive)}
#'      \item{...}{gene expression value in different samples}
#' }
#'
NULL
