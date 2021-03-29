#' @name validate_signature
#' @title Validate the signature model in a data sets
#' @description
#' @param
#' @import dplyr
#' @import purrr
#' @import cli
#' @export
#' @author
#' @return
#' @examples
#'
validate_signature <- function(train_sig = NULL, surv_data, id = "ids", time = "time", event = "status",
                               exp_data, labels_value = NULL, cut_p = NULL, keep_data = TRUE) {
  test_obj(obj = train_sig)

  ## get all signature genes
  all_beta <- purrr::reduce(train_sig$beta_value, bind_rows)

  ## check genes whether all gene in the exp_data
  if (!all(is.element(all_beta$Variable, rownames(exp_data)))) {
    last_genes <- setdiff(all_beta$Variable, rownames(exp_data))
    cli::cli_alert_info("These genes not in the expression profile: {last_genes}")
    stop("Not all genes of signatures in the expression profile!")
  }

  ## logrank analysis
  cli::cli_process_start("Checking the data")
  stopifnot(is.data.frame(surv_data))
  stopifnot(all(is.element(c(id, time, event), colnames(surv_data))))
  stopifnot(is.matrix(exp_data), is.numeric(exp_data))
  cli::cli_process_done()

  # get overlap samples
  cli::cli_process_start("Getting overlap samples")
  co_samples <- intersect(surv_data[[id]], colnames(exp_data))
  surv_data <- surv_data[match(co_samples, surv_data[[id]]) ,]
  exp_data <- exp_data[, co_samples]
  cli::cli_process_done()

  # get beta value
  cli::cli_process_start("Getting beta value")
  beta <- all_beta
  cli::cli_process_done()

  # get risk score
  cli::cli_process_start("Calculating risk score")
  signature_list <- purrr::map(train_sig$beta_value, dplyr::pull, Variable)
  all_score <- purrr::map(signature_list, function(x) {
    risk_score(exp_data = exp_data, genes = x, beta = beta) %>%
      dplyr::left_join(surv_data, by = id)
    }) %>%
    purrr::reduce(bind_rows) %>%
    dplyr::group_by(signature) %>%
    nest() %>%
    ungroup()
  cli::cli_process_done()

  # set labels
  cli::cli_process_start("Setting labels")
  case <- all_score %>%
    dplyr::pull(data) %>%
    purrr::map(~label_sample(score = .x$risk_score, value = labels_value))
  labeled_sample <- all_score %>%
    dplyr::mutate(data = purrr::map2(data, case, function(x, y) {x %>% dplyr::mutate(labels = y$labeled_sample)}),
                  cutoff_value = purrr::map_dbl(case, ~.x$value))
  cli::cli_process_done()

  # do logrank test
  cli::cli_process_start("Evaluating logrank test")
  logrank_res <- labeled_sample %>%
    dplyr::mutate(logrank_pval = purrr::map_dbl(data, function(x) {
      logrank_p(data = x, time = time, event = event, variate = "labels", verbose = F) %>% pull(p_value)
    }) %>% unlist()) %>%
    dplyr::select(signature, cutoff_value, logrank_pval, data)
  cli::cli_process_done()

  # get result
  cli::cli_process_start("Filtering result")
  if (is.null(cut_p)) {
    if (isTRUE(keep_data)) {
      res <- logrank_res
    } else {
      res <- dplyr::select(logrank_res, -data)
    }
  } else {
    if (isTRUE(keep_data)) {
      res <- dplyr::filter(logrank_res, logrank_pval < cut_p)
    } else {
      res <- dplyr::filter(logrank_res, logrank_pval < cut_p) %>%
        dplyr::select(-data)
    }
  }
  class(res) <- c("validate_signature", class(res))
  cli::cli_process_done()

  return(res)
}

#' @name validate_unicox
#' @title Do the univariate cox analysis in validation data sets.
#' @param obj the `validate_signature` object get from \link{validate_signature}() .
#' @param type Use which variate to do the univariate cox analysis, if `continuous`, use the `risk_score`;
#' if `discrete`, use the `labels`.
#' @param cut_p the cutoff p value of univariate cox analysis. Default 0.05.
#' @importFrom dplyr mutate filter select pull
#' @importFrom purrr map
#' @importFrom cli cli_process_start cli_process_done
#' @return return a `validate_signature` object with `unicox_pval` column.
#' @export
#' @example
#' \dontrun{
#'   uni_cox_obj <- validate_unicox(obj, type = "discrete", cut_p = 1)
#' }
#'

validate_unicox <- function(obj, type = c("continuous", "discrete"), cut_p = 0.05) {

  validate_obj(obj)
  cli::cli_process_start("Doing univariate cox analysis")
  var <- ifelse(match.arg(type) == "continuous", "risk_score", "labels")
  obj %<>% dplyr::mutate(unicox_pval = purrr::map(data, optimal_cox, variate = var, multicox = FALSE,
                                                  global_method = "wald")) %>%
    dplyr::mutate(unicox_pval = purrr::map(unicox_pval, dplyr::pull, p_value) %>% unlist()) %>%
    dplyr::filter(unicox_pval < cut_p)

  if (rlang::has_name(obj, "multicox_pval")) {
    res <- obj %>% select(1:3, 6, 4:5)
  } else {
    res <- obj %>% select(1:3, 5, 4)
  }
  cli::cli_process_done()

  return(res)
}


#' @name validate_multicox
#' @title Do the multivariate cox analysis in validation data sets.
#' @param obj the `validate_signature` object get from \link{validate_signature}() .
#' @param type Use which variate to do the multivariate cox analysis, if `continuous`, use the `risk_score`;
#' if `discrete`, use the `labels`.
#' @param cut_p the cutoff p value of multivariate cox analysis. Default 0.05.
#' @importFrom dplyr mutate filter select pull
#' @importFrom purrr map
#' @importFrom cli cli_process_start cli_process_done
#' @return return a `validate_signature` object with `multicox_pval` column.
#' @export
#' @examples
#' \dontrun{
#'   multi_cox_obj <- validate_multicox(obj = uni_cox_obj, type = "discrete", covariate = c("Age", "Gender"), cut_p = 1)
#'   multi_cox_obj <- validate_multicox(obj = obj, type = "discrete", covariate = c("Age", "Gender"), cut_p = 1)
#'   ## covariate = NULL
#'   multi_cox_obj <- validate_multicox(obj = uni_cox_obj, type = "discrete", cut_p = 1)
#' }

validate_multicox <- function(obj, type = c("continuous", "discrete"), covariate = NULL, cut_p = 0.05) {

  validate_obj(obj)
  ## type
  type <- match.arg(type)

  cli::cli_process_start("Doing Multivariate cox analysis")
  if (is.null(covariate)) {
    cli::cli_alert_info("The `covariate` is NULL, keep univariate result!")
    ## check whether done univariate cox analysis
    if (rlang::has_name(obj, "unicox_pval")) {
      res <- obj
    } else {
      res <- train_unicox(obj, type = type, cut_p = cut_p)
    }
  } else {
    stopifnot(is.element(covariate, colnames(obj$data[[1]])))
    uni_type <- ifelse(match.arg(type) == "continuous", "risk_score", "labels")
    covars <- c(uni_type, covariate)
    obj %<>% dplyr::mutate(multicox_pval = purrr::map(data, optimal_cox, variate = covars, multicox = TRUE,
                                                      global_method = "wald")) %>%
      dplyr::mutate(multicox_pval = purrr::map(multicox_pval, function(x) {
        x %>% filter(stringr::str_detect(Variable, uni_type)) %>% select(p_value)
      }) %>% unlist()) %>%
      dplyr::filter(multicox_pval < cut_p)
    if (rlang::has_name(obj, "unicox_pval")) {
      res <- obj %>% select(1:4, 6, 5)
    } else {
      res <- obj %>% select(1:3, 5, 4)
    }
  }
  cli::cli_process_done()

  return(res)
}



### useful function ######
validate_obj <- function(obj) {
  if (!is(obj, "validate_signature")) {
    stop("The `obj` is not a `training_signature` object!")
  } else {
    if (!rlang::has_name(obj, "data")) {
      stop("There is no information to do the analysis, please set the `keep_data` to TRUE in the
           `validate_signature` function to get the needed data!")
    }
  }
}
