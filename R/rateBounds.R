#' Get rate bounds
#'
#' Get vector of min and max rates
#' @param model rateModel object
#' @rdname getRateBounds
#' @name getRateBounds
#' @return vector of parameter values
#' @include rateModel-class.R
#' @export 
methods::setGeneric("getRateBounds", function(model) {
  standardGeneric("getRateBounds")
})

#' @name getRateBounds
#' @rdname getRateBounds
#' @aliases getRateBounds,rateModel,rateModel-method 
methods::setMethod("getRateBounds", signature(model = "rateModel"), function(model) {
  return(model@phylogeny$getRateBounds())
})

#' Set rate bounds
#'
#' Set rate bounds
#' @param model rateModel object
#' @param rMin rate lower bound
#' @param rMax rate upper bound
#' @rdname setRateBounds
#' @name setRateBounds
#' @return vector of parameter values
#' @include rateModel-class.R
#' @export 
methods::setGeneric("setRateBounds", function(model,rMin=NA,rMax=NA) {
  standardGeneric("setRateBounds")
})

#' @name setRateBounds
#' @rdname setRateBounds
#' @aliases setRateBounds,rateModel,rateModel-method 
methods::setMethod("setRateBounds", signature(model = "rateModel"), function(model,rMin=NA,rMax=NA) {
  rDefault=getRateBounds(model)
  rVec=c(rMin,rMax)
  ## Any NA elements get replaced with existing values
  rVec[is.na(rVec)]=rDefault[is.na(rVec)]
  ## If any elements are not finite throw an error
  if(!all(is.finite(rVec))){
    stop("All rate bounds must be finite numeric values")  
  }
  return(model@phylogeny$setRateBounds(rVec[1],rVec[2]))
})
