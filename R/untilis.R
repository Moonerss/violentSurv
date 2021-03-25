#' @name logrank_p
#' @title calculate logrank test p value
#' @description calculate logrank test p value
#' @param data a `data.frame` containing variables, time and status.
#' @param time column name specifying time, default is 'time'.
#' @param event column name specifying event status, default is 'status'.
#' @param variate column names of specific variables in the survival model.
#' @param verbose output other useful information, default TRUE.
#' @note If there is some wrong with the logrank test, the function will set logrank p value to 1 automatically.
#'   So set \code{verbose} argument to TRUE is necessary to check these events.
#' @importFrom survival survdiff
#' @importFrom stats pchisq as.formula
#' @importFrom tibble tibble
#' @import cli
#' @return return a tibble object. The first column is the variate name,
#' the second column is the logrank test p value, the third column is the Survival Object.
#' @export
#' @examples
#' library(survival)
#' data(ovarian)
#' logrank_p(data = ovarian, time = "futime", event = "fustat", variate = "rx")
#' logrank_p(data = ovarian, time = "futime", event = "fustat", variate = c("rx", "resid.ds"))
logrank_p <- function(data, time = "time", event = "status", variate, verbose = TRUE) {
  options(stringAsFactors = FALSE)

  fm <- paste0("Surv(", time, ", ", event, ") ~ ",
               ifelse(length(variate) == 1, variate, paste(variate, collapse = " + ")))

  survObject <- tryCatch(survival::survdiff(as.formula(fm), data = data),
                         error = function(e) {
                           if (verbose) {
                             cli::cli_alert_danger(paste("some problems with:", paste(variate, collapse = "; ")))
                             cli::cat_bullet("Set the p value to: 1")
                           }
                         })
  if (is.null(survObject)) {
    p.value <- 1
  } else {
    p.value <- stats::pchisq(survObject$chisq, length(survObject$n) - 1, lower.tail = F)
  }

  res <- tibble::tibble(Variable = paste(variate, collapse = "/"), p_value = p.value, models = as.character(fm))
  return(res)
}


#' open multiple threads process
#' @param nThreads The number of threads
#' @param maxSize The max memory size in global, default 500MB, the unit is MB
#' @param verbose output other useful information
#' @import future
#' @import cli
#' @importFrom future availableCores plan
#' @export
#' @examples
#' \dontrun{
#'   enableParallel(nThreads = 3)
#' }
enableParallel <- function(nThreads = NULL, maxSize = 500, verbose = FALSE) {
  nCores <- future::availableCores()
  options(future.globals.maxSize = maxSize*1024^2)
  if (verbose) cli::cli_alert_info(paste0("The maxSize is ", maxSize, " Mb."))
  if (is.null(nThreads)) {
    if (nCores < 4) {
      nThreads <- nCores
    } else {
      nThreads <- nCores - 2
    }
  }
  if (!is.numeric(nThreads) || nThreads < 2) {
    stop("nThreads must be numeric and at least 2.", call. = F)
  } else {
    if (nThreads > nCores) {
      if (verbose) {
        cli::cli_alert_info("Requested number of threads is higher than number of available processors (or cores)")
        cli::cli_alert_info(paste("The max working processes is", nCores))
      }
    }
    if (verbose) cli::cli_alert_info(paste("Allowing parallel execution with up to", nThreads, "working processes."))
    future::plan("multicore", workers = nThreads)
    invisible(nThreads)
  }
}
