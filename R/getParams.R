#' Get params
#'
#' Returns a vector of parameter values
#' @param x rateModel object
#' @rdname getParams
#' @name getParams
#' @return vector of parameter values
#' @include rateModel-class.R
#' @export 
methods::setGeneric("getParams", function(obj) {
  standardGeneric("getParams")
})

#' @name getParams
#' @rdname getParams
#' @aliases getParams,rateModel,rateModel-method 
methods::setMethod("getParams", signature(obj = "rateModel"), function(obj) {
  return(obj@phylogeny$getParams())
})