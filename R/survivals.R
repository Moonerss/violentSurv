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
#' @importFrom  dplyr bind_cols left_join group_by ungroup mutate filter select all_of
#' @importFrom  purrr map reduce map_dbl map2
#' @importFrom tidyr nest
#' @export
#' @return
#' data(blca_clinical)
#' data(blca_exp)
#' train_signature(surv_data = blca_clinical, id = "ids", time = "time", event = "status",
#'                 exp_data = blca_exp, genes = c("TP53", "BRCA1", "BRCA2", "RB1"), num = 2:4,
#'                 beta = NULL, labels_value = 1, cut_p = 1, keep_data = F)

train_signature <- function(surv_data, id = "ids", time = "time", event = "status",
                            exp_data, genes, num, beta = NULL, labels_value = NULL,
                            cut_p = NULL, keep_data = FALSE) {
  stopifnot(is.data.frame(surv_data))
  stopifnot(all(is.element(c(id, time, event), colnames(surv_data))))
  stopifnot(is.matrix(exp_data), is.numeric(exp_data), all(is.element(genes, rownames(exp_data))))

  # get overlap samples
  co_samples <- intersect(surv_data[[id]], colnames(exp_data))
  surv_data <- surv_data[match(co_samples, surv_data[[id]]) ,]
  exp_data <- exp_data[, co_samples]

  # get beta value
  if (is.null(beta)) {
    dat <- surv_data[, c(id, time, event)] %>%
      dplyr::bind_cols(as.data.frame(t(exp_data[genes,])))
    cox_obj <- run_cox(data = dat, time = time, event = event,
                       variate = genes, multicox = F, global_method = "wald")
    beta <- get_beta(cox_obj)
  } else {
    beta <- beta
  }

  # get signature combinations
  signature_list <- combn_signature(genes = genes, n = num)

  # get risk score
  all_score <- purrr::map(signature_list, function(x) {
    x %<>% t() %>% as.data.frame()
    purrr::map(x, function(y) {
      risk_score(exp_data = exp_data, genes = y, beta = beta) %>%
        dplyr::left_join(dplyr::select(surv_data, dplyr::all_of(id), dplyr::all_of(time), dplyr::all_of(event)), by = id)
      })
    }) %>%
    purrr::map(~purrr::reduce(.x, bind_rows)) %>%
    purrr::reduce(bind_rows) %>%
    dplyr::group_by(signature) %>%
    nest() %>%
    ungroup()

  # get labels
  case <- all_score %>%
    dplyr::pull(data) %>%
    purrr::map(~label_sample(score = .x$risk_score, value = labels_value))
  labeled_sample <- all_score %>%
    dplyr::mutate(data = purrr::map2(data, case, function(x, y) {x %>% dplyr::mutate(labels = y$labeled_sample)}),
                  cutoff_value = purrr::map_dbl(case, ~.x$value))

  # do logrank test
  logrank_res <- labeled_sample %>%
    dplyr::mutate(logrank_pval = purrr::map_dbl(data, function(x) {
      tryCatch(logrank_p(data = x, time = time, event = event, variate = "labels"),
               error = function(e) {2})
    }) %>% unlist()) %>%
    dplyr::select(signature, cutoff_value, logrank_pval, data)

  # get result
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
  return(res)
}
