#' Set params
#'
#' Set parameter values for any length of parameters
#' @param obj rateModel object
#' @param x parameter values
#' @param index an integer vector of parameter indicies to be changed (note: parameters are 0 indexed)
#' @rdname setParams
#' @name setParams
#' @include rateModel-class.R
#' @export 
methods::setGeneric("setParams", function(obj,x,index) {
  standardGeneric("setParams")
})

#' @name setParams
#' @rdname setParams
#' @aliases setParams,rateModel,rateModel-method 
methods::setMethod("setParams", signature(obj = "rateModel"), function(obj,x,index) {
  return(obj@phylogeny$setParams(x,index))
})