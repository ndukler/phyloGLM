#' Compute parameter gradients
#'
#' Compute parameter gradients for given rate model object
#' @param x a set of parameters around which to compute the gradient (optional)
#' @param model rateModel object
#' @param scale value to scale log-likelihood by (default=1)
#' @param threads number of threads to use
#' @return vector of gradients
#' @name phyloGrad
#' @include rateModel-class.R
#' @rdname phyloGrad
methods::setGeneric("phyloGrad", function(x, model, scale = 1, eps = .Machine$double.eps, threads = 1, index = NULL) {
  standardGeneric("phyloGrad")
})


#' @name phyloGrad
#' @rdname phyloGrad
methods::setMethod("phyloGrad", signature(x = "missing", model = "rateModel"), function(x, model, scale = 1,
                                                                                        eps = .Machine$double.eps, threads = 1, index = NULL) {
  g <- model@phylogeny$grad(
    model@alleleData$alleleData@data@x,
    model@rateDM@x,
    model@piDM@x, scale, threads
  )[!model@fixed]
  return(g)
})

#' @name phyloGrad
#' @rdname phyloGrad
methods::setMethod("phyloGrad", signature(x = "numeric", model = "rateModel"), function(x, model, scale = 1,
                                                                                        eps = .Machine$double.eps, threads = 1, index = NULL) {
  ## message(paste("X:",paste(x,collapse = "\t")))
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
  ## Set parameters
  setParams(model, x, index)
  ## Evaluate gradient
  g <- phyloGrad(model = model, scale = scale, eps = eps, threads = threads)
  ## Restore original parameter values
  setParams(model, initP, 0:(length(initP) - 1))
  ## message(paste("G:",paste(g,collapse = "\t")))
  return(g)
})
