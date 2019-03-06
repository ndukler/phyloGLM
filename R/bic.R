#' Compute BIC for rate models
#'
#' Computes the bayesian information criterion (BIC) for a given rateModel.
#' Note that when comparing two models, the one with the smaller BIC is preferred. 
#' @param model rateModel object
#' @name bic
#' @include rateModel-class.R
#' @rdname bic
#' @return numeric BIC value
#' @examples
#' 
#' @export
methods::setGeneric("bic", function(model,i) {
  standardGeneric("bic")
})

#' @name bic
#' @rdname bic
methods::setMethod("bic", signature(model = "rateModel"), function(model) {
  l=ll(model=model)
  k=sum(!model@fixed)
  nSites=model@rateDM@nrow
  return((log(nSites)*k)-(2*l))
})