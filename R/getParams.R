#' Plot coefficents
#'
#' Returns a vector of parameter values
#' @param model rateModel object
#' @rdname getParams
#' @name getParams
#' @return vector of parameter values
#' @include rateModel-class.R
#' @export
methods::setGeneric("getParams", function(model) {
  standardGeneric("getParams")
})

#' @name getParams
#' @rdname getParams
#' @aliases getParams,rateModel,rateModel-method
methods::setMethod("getParams", signature(model = "rateModel"), function(model) {
  return(model@phylogeny$getParams())
})
