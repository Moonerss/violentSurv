#' @name combn_signature
#' @title get all combination of specific number
#' @description get all combination of specific number
#' @param genes the candidate set of genes
#' @param n the number of genes to select
#' @importFrom  utils combn
#' @importFrom  purrr map
#' @export
#' @return return a list of combination, each row is a combination
#' @examples
#' combn_signature(letters[1:4], n = 2)
#' combn_signature(letters[1:4], n = 2:3)
#'
combn_signature <- function(genes, n) {
  res <- purrr::map(n, ~t(combn(genes, m = .x)))
  return(res)
}

#' @name get_beta
#' @title beta coefficient of univariate cox analysis
#' @description get the beta coefficient of univariate cox analysis
#' @param run_cox_obj the `run_cox` object get from the \link{run_cox}()
#' @importFrom dplyr select
#' @export
#' @return return a `data.frame`, the first column is `Variable`, the second column is `beta`.
#' @examples
#' library(survival)
#' data(lung)
#' cox_obj <- run_cox(data = lung, time = "time", event = "status",
#'   variate = c("age", "sex", "ph.ecog"), multicox = FALSE, global_method = "wald")
#' get_beta(cox_obj)
#'
get_beta <- function(run_cox_obj) {
  run_cox_obj %>% dplyr::select(Variable, beta)
}


#' @name risk_score
#' @title calculate the expression risk score
#' @description batch calculate the expression risk score
#' @details this function calculate the risk score based on expression profile.
#' `risk score = gene_1 * coef_1 + ... + gene_n * coef_n`
#' @param exp_data a `matrix` of expression profile.
#' @param genes a character vector contain the genes in the signature.
#' @param beta a  contained the beta value
#' correspond to the genes.
#' @importFrom  dplyr pull filter tibble
#' @export
#' @return return a `tibble`
#' @examples
#' \dontrun{
#'     data(meta_dat)
#'     data(blca_exp)
#'     genes <- colnames(meta_dat)[4:10]
#'     cox_res <- run_cox(data = meta_dat, time = "time", event = "status",
#'                        variate = genes, multicox = FALSE, global_method = "wald")
#'     beta <- get_beta(cox_res)
#'     scores <- risk_score(exp_data = blca_exp, genes, beta)
#' }
#'
risk_score <- function(exp_data, genes, beta) {
  stopifnot(is.matrix(exp_data), is.numeric(exp_data))
  stopifnot(is.vector(genes), is.character(genes))
  stopifnot(all(is.element(genes, dplyr::pull(beta, Variable))))
  stopifnot(all(is.element(genes, rownames(exp_data))))

  # extract the genes expression
  exp_data <- exp_data[genes,]

  # extract the gene beta
  sub_beta <- beta %>% dplyr::filter(is.element(Variable, genes)) %>%
    {.[match(genes, dplyr::pull(.,Variable)),]}


  # weighted expression
  weighted_exp <- exp_data * dplyr::pull(sub_beta, beta)
  scores <- apply(weighted_exp, 2, sum)

  res <- dplyr::tibble(
    signature = paste(genes, collapse = " "),
    # beta_value = tibble::as_tibble(sub_beta),
    ids = names(scores),
    risk_score = scores
  )

  return(res)
}


#' @name label_sample
#' @title High and Low
#' @description Differentiate samples based on risk score
#' @details If the score is higher than the specific value, the samples are labeled
#' as `High_risk`, otherwise as `Low_risk`
#' @param score the vector of risk score.
#' @param value the specific value to distinguish different samples. Default is `NULL`,
#' we used the median of risk score to divide samples.
#' @importFrom stats median
#' @export
#' @examples
#' label_sample(1:5)
#' label_sample(1:5, value = 3)
label_sample <- function(score, value = NULL) {
  if (is.null(value)) {
    labeled_sample <- ifelse(score > median(score), "High_risk", "Low_risk")
    value <- median(score)
  } else {
    labeled_sample <- ifelse(score > value, "High_risk", "Low_risk")
    value <- value
  }
  return(list(labeled_sample = labeled_sample, value = value))
}

