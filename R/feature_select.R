#' @description: learn models for elastic net, lasso, ridge and adaptive elastic net
#' @title learn_glmnet
#' @param gene_data matrix of expression data (row: gene, column: patients)
#' @param y_cox Surv object
#' @param method "elastic" for elastic net, "lasso" for lasso, "ridge" for ridge and "adaptive" for adaptive elastic net
#' @param nfold number of nfold for cross-validation - default is 10
#' @param lambda "lambda.min" or "lambda.1se"
#'
#' @export
#'
learn_glmnet <- function(gene_data, y_cox, method = c("lasso", "elastic", "ridge", "adaptive"),
                         nfold = 10, alpha = 0.3, lambda = c("lambda.min", "lambda.1se")) {

  ## set the result
  cv_fits <- list()

  ## get options
  type <- match.arg(type)
  lambda <- match.arg(lambda)

  ## run feature select
  if (method == "lasso") {
    message("Selecting features by lasso ...")
    cvfit <- cv.glmnet(gene_data, y_cox, family = "cox", alpha = 1, nfolds = nfold, grouped = T, standardize = F)
  } else if (method == "elastic") {
    message("Selecting features by elastic net ...")
    cvfit <- cv.glmnet(gene_data, y_cox, family = "cox", alpha = alpha, nfolds = nfold, grouped = T, standardize = F)
  } else if (method == "ridge") {
    message("Selecting features by ridge net ...")
    cvfit <- cv.glmnet(gene_data, y_cox, family = "cox", alpha = 0, nfolds = nfold, grouped = T, standardize = F)
  } else if (method == "adaptive") {
    message("Selecting features by adaptive elastic net ...")
    fit_ridge <- cv.glmnet(gene_data, y_cox, family = "cox", alpha = 0, nfolds = nfold, grouped = T, standardize = F)
    betas_ridge <- as.numeric(coef(fit_ridge, "lambda.min"))
    w_AEN <- abs(1 / betas_ridge)
    cvfit <- cv.glmnet(gene_data, y_cox, family = "cox", nfolds = nfold, standardize = F, alpha = alpha,
                       grouped = T, penalty.factor = w_AEN)
  }

  ## get selected featires
  myCoefs <- coef(cvfit, s = lambda)
  cv_fits$res <- tibble::tibble(
    features = row.names(myCoefs)[myCoefs@i[-1]],
    coef = myCoefs@x[myCoefs@i[-1]]
  )
  cv_fits$cvfit <- cvfit

  return(cv_fits)
}
