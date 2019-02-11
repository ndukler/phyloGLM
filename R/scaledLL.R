#' Compute scaled negative log-likelihood
#'
#' Compute scaled log-likelihood for given rate model object
#' @param x a set of parameters to compute the log-likelihood for (optional)
#' @param obj rateModel object
#' @param scale value to scale log-likelihood by (default=1)
#' @param threads number of threads to use
#' @return scaled log-likelihood
#' @name scaledLL
#' @include rateModel-class.R
#' @rdname scaledLL
methods::setGeneric("scaledLL", function(x,obj,scale=1,threads=1,...) {
  standardGeneric("scaledLL")
})


#' @name scaledLL
#' @rdname scaledLL
methods::setMethod("scaledLL", signature(x="missing",obj = "rateModel"), function(x,obj,scale=1,threads=1) {
  return(sum(siteLL(obj,threads=threads))*scale)
})

#' @name scaledLL
#' @rdname scaledLL
methods::setMethod("scaledLL", signature(x="numeric",obj = "rateModel"), function(x,obj,scale=1,threads=1) {
  setParams(obj,x,which(!obj@fixed)-1)
  return(sum(siteLL(obj,threads = threads))*scale)
})
