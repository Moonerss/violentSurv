#' @name logrank_p
#' @title calculate logrank test p value
#' @description calculate logrank test p value
#' @param data a `data.frame` containing variables, time and status.
#' @param time column name specifying time, default is 'time'.
#' @param event column name specifying event status, default is 'status'.
#' @param variate column names of specific variables in the survival model.
#' @importFrom survival survdiff
#' @importFrom stats pchisq as.formula
#' @return return the logrank test p value
#' @export
#' @examples
#' library(survival)
#' data(ovarian)
#' logrank_p(data = ovarian, time = "futime", event = "fustat", variate = "rx")
#' logrank_p(data = ovarian, time = "futime", event = "fustat", variate = c("rx", "resid.ds"))
logrank_p <- function(data, time = "time", event = "status", variate) {
  options(stringAsFactors = FALSE)

  data$time <- data[[time]]
  data$event <- data[[event]]

  fm <- as.formula(ifelse(length(variate) == 1,
               paste("Surv(time, event)~", variate),
               paste("Surv(time, event)~", paste(variate, collapse = "+"))))

  survObject <- tryCatch(survival::survdiff(fm, data = data),
                         error = function(e) {
                           stop(paste("some problems with:", paste(variate, collapse = "; ")))
                         })
  p.value <- stats::pchisq(survObject$chisq, length(survObject$n) - 1, lower.tail = F)

  return(p.value)
}

