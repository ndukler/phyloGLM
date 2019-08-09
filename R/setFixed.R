#' Set fixed
#'
#' Set whether parmeters aree fixed
#' @param obj rateModel object
#' @param index a 0-based integer vector of parameter indicies to be changed (note: parameters are 0 indexed)
#' @param x logical vector of whether parameter values are fixed
#' @rdname setFixed
#' @name setFixed
#' @include rateModel-class.R

methods::setGeneric("setFixed", function(obj, x, index) {
  standardGeneric("setFixed")
})

#' @name setParams
#' @rdname setParams
#' @aliases setParams,rateModel,rateModel-method
methods::setMethod("setFixed", signature(obj = "rateModel"), function(obj, x, index) {
  updateFixed(obj@fixed, x, index)
})
