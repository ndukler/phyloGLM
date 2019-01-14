#' Compute scaled negative log-likelihood
#'
#' Compute scaled log-likelihood for given rate model object
#' @param x a set of parameters to compute the log-likelihood for (optional)
#' @param obj rateModel object
#' @return scaled log-likelihood
#' @name scaledLL
#' @include rateModel-class.R
#' @rdname scaledLL
methods::setGeneric("scaledLL", function(x,obj,scale,...) {
  standardGeneric("scaledLL")
})


#' @name scaledLL
#' @rdname scaledLL
methods::setMethod("scaledLL", signature(x="missing",obj = "rateModel"), function(x,obj,scale) {
  
  return(ll)
})

#' @name scaledLL
#' @rdname scaledLL
methods::setMethod("scaledLL", signature(x="numeric",obj = "rateModel"), function(x,obj,scale) {
  ll=logLikelihood(obj = obj)
  return(ll)
})
