#' Compute hessian
#'
#' Compute parameter hessian for given rate model object
#' @param x a set of parameters around which to compute the gradient (optional)
#' @param model rateModel object
#' @param scale value to scale log-likelihood by (default=1)
#' @param threads number of threads to use
#' @return vector of gradients
#' @name phyloHess
#' @include rateModel-class.R
#' @rdname phyloHess
methods::setGeneric("phyloHess", function(x, model, scale = 1, eps = .Machine$double.eps, threads = 1) {
  standardGeneric("phyloHess")
})


#' @name phyloHess
#' @rdname phyloHess
methods::setMethod("phyloHess", signature(x = "missing", model = "rateModel"), function(x, model, scale = 1,
                                                                                        eps = .Machine$double.eps, threads = 1) {
  ## Get initial parameters
  x <- getParams(model)
  ## Calculate stepsize for optimal finite gradient approximation
  h <- abs((eps^(1 / 3)) * x)
  ## Set lower bound on stepsize at 10^-14
  h[h < 10^-14] <- 10^-14
  ## Initialize likelihood matrices
  f_plus <- matrix(nrow = length(x), ncol = length(x))
  f_minus <- matrix(nrow = length(x), ncol = length(x))
  ## Initialize vector for realized h
  h_real <- matrix(nrow = length(x), ncol = length(x))
  for (p in 1:length(x)) {
    ## Compute alternate values parameters
    temp_plus <- x[p] + h[p]
    temp_minus <- x[p] - h[p]
    ## Compute realized h
    h_real[, p] <- temp_plus - temp_minus
    ## Do plus calculation
    setParams(model, temp_plus, p - 1)
    f_plus[, p] <- phyloGLM:::phyloGrad(model = model, scale = scale, eps = eps, threads = threads)
    ## Minus calculation
    setParams(model, temp_minus, p - 1)
    f_minus[, p] <- phyloGLM:::phyloGrad(model = model, scale = scale, eps = eps, threads = threads)
    ## Reset parameter
    setParams(model, x[p], p - 1)
  }
  ## Compute Gradients
  h <- (f_plus - f_minus) / h_real
  return(h)
})

#' @name phyloHess
#' @rdname phyloHess
methods::setMethod("phyloHess", signature(x = "numeric", model = "rateModel"), function(x, model, scale = 1,
                                                                                        eps = .Machine$double.eps, threads = 1) {
  setParams(model, x, which(!model@fixed) - 1)
  return(phyloHess(model = model, scale = scale, eps = eps, threads = threads))
})
