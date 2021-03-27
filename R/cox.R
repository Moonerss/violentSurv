#' @name run_cox
#' @title Run Cox Analysis in Batch Mode
#' @param data a `data.frame` containing variables, time and event.
#' @param time column name specifying time, default is 'time'.
#' @param event column name specifying event status, default is 'status'.
#' @param variate column names specifying variables.
#' @param multicox default `FALSE`. If `TRUE`, do the multiple variate cox analysis.
#' @param global_method method used to obtain global p value for cox model,
#' should be one of "likelihood", "wald", "logrank".
#' The likelihood-ratio test, Wald test, and score logrank statistics.
#' These three methods are asymptotically equivalent. Default is "wald".
#' @importFrom  survival coxph
#' @importFrom  dplyr tibble
#' @importFrom  purrr map reduce
#' @importFrom  magrittr %<>%
#' @return return a `run_cox` object including the summary information of cox analysis
#' @export
#' @examples
#' library(survival)
#' data(lung)
#' run_cox(data = lung, time = "time", event = "status",
#'   variate = c("age", "sex", "ph.ecog"), multicox = FALSE, global_method = "wald")
#' run_cox(data = lung, time = "time", event = "status",
#'   variate = c("age", "sex", "ph.ecog"), multicox = TRUE, global_method = "wald")
#'
run_cox <- function(data, time = "time", event = "status",
                    variate, multicox = FALSE,
                    global_method = c("likelihood", "wald", "logrank")) {
  stopifnot(is.data.frame(data))
  stopifnot(is.character(variate), is.vector(variate), all(is.element(variate, colnames(data))))

  ## choose univariate cox or multiple variate cox analysis
  if (isTRUE(multicox)) {
    if (length(variate) < 2) {
      stop("For multiple variate cox analysis, there must be more than one variate!")
    } else {
      res <- one_cox(var = variate, data, time, event, multicox, method = global_method)
    }
  } else {
    res <- purrr::map(variate, one_cox, data, time, event, multicox, method = global_method) %>%
      purrr::reduce(dplyr::add_row)
  }

  class(res) <- c("run_cox", class(res))
  return(res)
}


#' @name run_cox_parallel
#' @title Run Cox Analysis in Batch Mode
#' @param data a `data.frame` containing variables, time and event.
#' @param time column name specifying time, default is 'time'.
#' @param event column name specifying event status, default is 'status'.
#' @param variate column names specifying variables.
#' @param multicox default `FALSE`. If `TRUE`, do the multiple variate cox analysis.
#' @param global_method method used to obtain global p value for cox model,
#' should be one of "likelihood", "wald", "logrank".
#' The likelihood-ratio test, Wald test, and score logrank statistics.
#' These three methods are asymptotically equivalent. Default is "wald".
#' @param verbose output other useful information.
#' @importFrom  survival coxph Surv
#' @importFrom  dplyr tibble
#' @importFrom  purrr reduce
#' @importFrom  magrittr %<>%
#' @import      furrr
#' @return return a `run_cox` object including the summary information of cox analysis
#' @export
#' @examples
#' library(survival)
#' data(lung)
#' run_cox_parallel(data = lung, time = "time", event = "status",
#'   variate = c("age", "sex", "ph.ecog"), multicox = FALSE, global_method = "wald")
#' run_cox_parallel(data = lung, time = "time", event = "status",
#'   variate = c("age", "sex", "ph.ecog"), multicox = TRUE, global_method = "wald")
#'
run_cox_parallel <- function(data, time = "time", event = "status",
                    variate, multicox = FALSE,
                    global_method = c("likelihood", "wald", "logrank"), verbose = TRUE) {
  stopifnot(is.data.frame(data))
  stopifnot(is.character(variate), is.vector(variate), all(is.element(variate, colnames(data))))

  ## choose univariate cox or multiple variate cox analysis
  if (isTRUE(multicox)) {
    if (length(variate) < 2) {
      stop("For multiple variate cox analysis, there must be more than one variate!", call. = FALSE)
    } else {
      res <- one_cox(var = variate, data, time, event, multicox, method = global_method)
    }
  } else {
    oplan <- future::plan()
    future::plan("multiprocess")
    on.exit(future::plan(oplan), add = TRUE)
    res <- furrr::future_map(variate, one_cox, data, time, event, multicox, method = global_method) %>%
      purrr::reduce(dplyr::add_row)
  }

  class(res) <- c("run_cox", class(res))
  return(res)
}

# one cox analysis
one_cox <- function(var, data, time = "time", event = "status", multicox = FALSE,
                    method = c("likelihood", "wald", "logrank")) {

  fm_chara <- paste0("Surv(", time, ", ", event, ") ~ ",
                     ifelse(multicox, paste(var, collapse = " + "), var))
  fm <- as.formula(fm_chara)

  cox <- tryCatch(coxph(fm, data = data),
                  error = function(e) {
                    stop("Something wrong with variable ", var)
                  })
  cox_res <- summary(cox)
  test_method <- switch(match.arg(method),
                        likelihood = "logtest",
                        wald =  "waldtest",
                        logrank = "sctest")
  res <- dplyr::tibble(
    Variable = rownames(cox_res$coefficients),
    Mulicox = multicox,
    beta = cox_res$coefficients[, "coef"],
    HR = cox_res$coefficients[, "exp(coef)"],
    lower_95 = cox_res$conf.int[, "lower .95"],
    upper_95 = cox_res$conf.int[, "upper .95"],
    p_value = cox_res$coefficients[, "Pr(>|z|)"],
    global_pval = cox_res[[test_method]]["pvalue"],
    models = fm_chara
  )

  return(res)
}
