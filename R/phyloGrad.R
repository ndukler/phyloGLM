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
  # ## Get initial parameters
  # x=getParams(model)
  # ## Calculate stepsize for optimal finite gradient approximation
  # h=abs((eps^(1/3))*x)
  # ## Set lower bound on stepsize at 10^-14
  # h[h<10^-14]=10^-14
  # ## Initialize likelihood vectors
  # f_plus=numeric(length(x))
  # f_minus=numeric(length(x))
  # ## Initialize vector for realized h
  # h_real=numeric(length(h))
  # for(p in 1:length(x)){
  #   ## Compute alternate values parameters
  #   temp_plus=x[p]+h[p]
  #   temp_minus=x[p]-h[p]
  #   ## Compute realized h
  #   h_real[p]=temp_plus-temp_minus
  #   ## Do plus calculation
  #   setParams(model,temp_plus,p-1)
  #   f_plus[p]=model@phylogeny$ll(model@alleleData$alleleData@data@x,model@rateDM@x,model@piDM@x,scale,threads)
  #   ## Minus calculation
  #   setParams(model,temp_minus,p-1)
  #   f_minus[p]=model@phylogeny$ll(model@alleleData$alleleData@data@x,model@rateDM@x,model@piDM@x,scale,threads)
  #   ## Reset parameter
  #   setParams(model,x[p],p-1)
  # }
  # ## Compute Gradients
  # g=(f_plus-f_minus)/h_reals
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
  return(g)
})
