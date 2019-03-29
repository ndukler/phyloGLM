#' Compute LRT for rate models
#'
#' Performs a likelihood ratio test for two models
#' @param h0 rateModel object that represents the null hypothesis
#' @param hA rateModel object that represents the alternative hypothesis (more complex model)
#' @name lrt
#' @include rateModel-class.R
#' @rdname lrt
#' @return a data table including D-statistic, the p-value, and the degrees of freedom
#' @examples
#' 
#' @export
methods::setGeneric("lrt", function(h0, hA) {
  standardGeneric("lrt")
})

#' @name lrt
#' @rdname lrt
methods::setMethod("lrt", signature(h0 = "rateModel", hA = "rateModel"), function(h0, hA) {
  nll_0 <- -ll(model = h0)
  nll_A <- -ll(model = hA)
  d <- 2 * (nll_0 - nll_A)
  deltaDF <- sum(!hA@fixed) - sum(!h0@fixed)
  p <- pchisq(d, lower.tail = FALSE, df = deltaDF)
  return(data.table::data.table(D = d, p = p, df = deltaDF))
})
