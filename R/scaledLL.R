#' Compute scaled negative log-likelihood
#'
#' Compute scaled log-likelihood for given rate model object
#' @param x a set of parameters to compute the log-likelihood for (optional). Only unfixed parameters.
#' @param model rateModel object
#' @param scale value to scale log-likelihood by (default=1)
#' @param threads number of threads to use
#' @param index integer vector of parameter indicies in vector (0-based)
#' @return scaled log-likelihood
#' @name scaledLL
#' @include rateModel-class.R
#' @rdname scaledLL
methods::setGeneric("scaledLL", function(x, model, scale = 1, threads = 1, index = NULL) {
  standardGeneric("scaledLL")
})


#' @name scaledLL
#' @rdname scaledLL
methods::setMethod("scaledLL", signature(x = "missing", model = "rateModel"), function(x, model, scale = 1, threads = 1, index = NULL) {
  if(!is.null(index)){
    warning("Index is specified without x and will be ignored.")
  }
  l <- model@phylogeny$ll(
    model@alleleData$alleleData@data@x,
    model@rateDM@x,
    model@piDM@x, scale, threads
  )
  if (is.nan(l)) {
    warning("NaN value for LL with parameters: ", paste(getParams(model), collapse = ","), "\n")
  } else if (is.infinite(l)) {
    warning("Inf value for LL with parameters: ", paste(getParams(model), collapse = ","), "\n")
  }
  return(l)
})

#' @name scaledLL
#' @rdname scaledLL
methods::setMethod("scaledLL", signature(x = "numeric", model = "rateModel"), function(x, model, scale = 1, threads = 1, index = NULL) {
  ## Save initial parameter values
  initP <- getParams(model)
  ## Default to index being all parameters, otherwise must specify index
  if(is.null(index)){
    if(length(x) != length(initP)){
      stop("If no index specified x must be same length as the full parameter vector")
    } else{
      index = 0:(length(initP)-1) 
    }
  } else {
    if(length(x) != length(index)){
      stop("length(x) does not equal length(index)")
    }
  }
  ## Set parameter values
  setParams(model, x, index)
  ## Evaluate LL
  l <- scaledLL(model = model, scale = scale, threads = threads)
  ## Restore original parameter values
  setParams(model, initP, 0:(length(initP) - 1))
  return(l)
})

#' Compute log-likelihood
#'
#' Compute the log-likelihood for given rate model object.
#' @param x a set of parameters to compute the log-likelihood for (optional)
#' @param model rateModel object
#' @param threads number of threads to use
#' @param index integer vector of parameter indicies in vector (0-based)
#' @return scaled log-likelihood
#' @name ll
#' @include rateModel-class.R
#' @rdname ll
#' @export
methods::setGeneric("ll", function(x, model, threads = 1, index = NULL) {
  standardGeneric("ll")
})


#' @name ll
#' @rdname ll
methods::setMethod("ll", signature(x = "missing", model = "rateModel"), function(x, model, threads = 1, index = NULL) {
  return(scaledLL(model = model, scale = 1, threads = threads, index = index))
})

#' @name ll
#' @rdname ll
methods::setMethod("ll", signature(x = "numeric", model = "rateModel"), function(x, model, threads = 1, index = NULL) {
  return(scaledLL(x = x, model = model, scale = 1, threads = threads, index = index))
})
