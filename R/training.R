#' @name train_signature
#' @title Train the signature model in a data sets
#' @description this function will filter the signature are not significant in the logrank test
#' @param surv_data a `data.frame` containing variables, ids, time and event.
#' @param id column name specifying samples ID , default is 'ids'.
#' @param time column name specifying time, default is 'time'.
#' @param event column name specifying event status, default is 'status'.
#' @param exp_data a `matrix` of expression profile. row names are genes, column names are sample ID.
#' @param genes a character vector of genes to construct risk model. They must be in `exp_data` row names.
#' @param num a numeric vector of the number of genes in the risk model.
#' @param beta if `NULL`, the function compute the beta value for each gene in univariate cox analysis;
#' or a data like the return object of \link{get_beta}()
#' @param labels_value the same like `value` in \link{label_sample}()
#' @param cut_p a p value to filter the result of `train_signature`, we only keep the risk model with p value
#' samller than the value
#' @param keep_data whether to keep the data used to compute logrank test in this function,
#' It may be very important for follow-up analysis.
#' @importFrom dplyr bind_cols left_join group_by ungroup mutate filter select all_of bind_rows
#' @importFrom purrr map reduce map_dbl map2
#' @importFrom tidyr nest
#' @import cli
#' @export
#' @examples
#' \dontrun{
#' data(blca_clinical)
#' data(blca_exp)
#' obj <- train_signature(surv_data = blca_clinical, id = "ids", time = "time", event = "status",
#'                        exp_data = blca_exp, genes = c("IL6", "TGFB3", "VHL", "CXCR4"), num = 2:4,
#'                        beta = NULL, labels_value = NULL, cut_p = 1, keep_data = T)
#' }
#'

train_signature <- function(surv_data, id = "ids", time = "time", event = "status",
                            exp_data, genes, num, beta = NULL, labels_value = NULL,
                            cut_p = NULL, keep_data = TRUE) {
  cli::cli_process_start("Checking the data")
  stopifnot(is.data.frame(surv_data))
  stopifnot(all(is.element(c(id, time, event), colnames(surv_data))))
  stopifnot(is.matrix(exp_data), is.numeric(exp_data), all(is.element(genes, rownames(exp_data))))
  cli::cli_process_done()

  # get overlap samples
  cli::cli_process_start("Getting overlap samples")
  co_samples <- intersect(surv_data[[id]], colnames(exp_data))
  surv_data <- surv_data[match(co_samples, surv_data[[id]]) ,]
  exp_data <- exp_data[, co_samples]
  cli::cli_process_done()

  # get beta value
  cli::cli_process_start("Getting beta value")
  if (is.null(beta)) {
    dat <- surv_data[, c(id, time, event)] %>%
      dplyr::bind_cols(as.data.frame(t(exp_data[genes,])))
    ## whether use parallel
    cox_obj <- optimal_cox(data = dat, time = time, event = event,
                           variate = genes, multicox = F, global_method = "wald")
    beta <- get_beta(cox_obj)
  } else {
    beta <- beta
  }
  cli::cli_process_done()

  # get signature combinations
  cli::cli_process_start("Combining signatures")
  signature_list <- combn_signature(genes = genes, n = num)
  cli::cli_process_done()

  # get risk score
  cli::cli_process_start("Calculating risk score")
  all_score <- purrr::map(signature_list, function(x) {
    x %<>% t() %>% as.data.frame()
    purrr::map(x, function(y) {
      risk_score(exp_data = exp_data, genes = y, beta = beta) %>%
        dplyr::left_join(surv_data, by = id)
      })
    }) %>%
    purrr::map(~purrr::reduce(.x, bind_rows)) %>%
    purrr::reduce(bind_rows) %>%
    dplyr::group_by(signature) %>%
    nest() %>%
    ungroup()
  all_score %<>% dplyr::mutate(beta_value = purrr::map(signature, function(x) {
    sig <- unlist(strsplit(x, split = " "))
    res <- dplyr::slice(beta, match(sig, Variable))
  })) %>% dplyr::mutate(beta_value = purrr::map(beta_value, function(x) {
    class(x) <- setdiff(class(x), "run_cox"); x}))
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
    dplyr::select(signature, cutoff_value, logrank_pval, beta_value, data)
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
  class(res) <- c("training_signature", class(res))
  cli::cli_process_done()

  return(res)
}


#' @name train_unicox
#' @title Do the univariate cox analysis in training data sets.
#' @param obj the `training_signature` object get from \link{train_signature}() .
#' @param type Use which variate to do the univariate cox analysis, if `continuous`, use the `risk_score`;
#' if `discrete`, use the `labels`.
#' @param cut_p the cutoff p value of univariate cox analysis. Default 0.05.
#' @importFrom dplyr mutate filter select pull
#' @importFrom purrr map
#' @importFrom cli cli_process_start cli_process_done
#' @return return a `training_signature` object with `unicox_pval` column.
#' @export
#' @example
#' \dontrun{
#'   uni_cox_obj <- train_unicox(obj, type = "discrete", cut_p = 1)
#' }
#'
train_unicox <- function(obj, type = c("continuous", "discrete"), cut_p = 0.05) {

  test_obj(obj)
  cli::cli_process_start("Doing univariate cox analysis")
  var <- ifelse(match.arg(type) == "continuous", "risk_score", "labels")
  obj %<>% dplyr::mutate(unicox_pval = purrr::map(data, optimal_cox, variate = var, multicox = FALSE,
                         global_method = "wald")) %>%
    dplyr::mutate(unicox_pval = purrr::map(unicox_pval, dplyr::pull, p_value) %>% unlist()) %>%
    dplyr::filter(unicox_pval < cut_p)

  if (rlang::has_name(obj, "multicox_pval")) {
    res <- obj %>% select(1:3, 7, 4, 5:6)
  } else {
    res <- obj %>% select(1:3, 6, 4:5)
  }
  cli::cli_process_done()

  return(res)
}


#' @name train_multicox
#' @title Do the multivariate cox analysis in training data sets.
#' @param obj the `training_signature` object get from \link{train_signature}() .
#' @param type Use which variate to do the multivariate cox analysis, if `continuous`, use the `risk_score`;
#' if `discrete`, use the `labels`.
#' @param cut_p the cutoff p value of multivariate cox analysis. Default 0.05.
#' @importFrom dplyr mutate filter select pull
#' @importFrom purrr map
#' @importFrom cli cli_process_start cli_process_done
#' @return return a `training_signature` object with `multicox_pval` column.
#' @export
#' @examples
#' \dontrun{
#'   multi_cox_obj <- train_multicox(obj = uni_cox_obj, type = "discrete", covariate = c("Age", "Gender"), cut_p = 1)
#'   multi_cox_obj <- train_multicox(obj = obj, type = "discrete", covariate = c("Age", "Gender"), cut_p = 1)
#'   ## covariate = NULL
#'   multi_cox_obj <- train_multicox(obj = uni_cox_obj, type = "discrete", cut_p = 1)
#' }

train_multicox <- function(obj, type = c("continuous", "discrete"), covariate = NULL, cut_p = 0.05) {

  test_obj(obj)
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
      res <- obj %>% select(1:4, 7, 5:6)
    } else {
      res <- obj %>% select(1:3, 6, 4:5)
    }
  }
  cli::cli_process_done()

  return(res)
}


### useful function ######
test_obj <- function(obj) {
  if (!is(obj, "training_signature")) {
    stop("The `obj` is not a `training_signature` object!")
  } else {
    if (!rlang::has_name(obj, "data")) {
      stop("There is no information to do the analysis, please set the `keep_data` to TRUE in the
           `train_signature` function to get the needed data!")
    }
  }
}
