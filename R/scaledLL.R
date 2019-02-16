#' Compute scaled negative log-likelihood
#'
#' Compute scaled log-likelihood for given rate model object
#' @param x a set of parameters to compute the log-likelihood for (optional)
#' @param model rateModel object
#' @param scale value to scale log-likelihood by (default=1)
#' @param threads number of threads to use
#' @return scaled log-likelihood
#' @name scaledLL
#' @include rateModel-class.R
#' @rdname scaledLL
methods::setGeneric("scaledLL", function(x,model,scale=1,threads=1,...) {
  standardGeneric("scaledLL")
})


#' @name scaledLL
#' @rdname scaledLL
methods::setMethod("scaledLL", signature(x="missing",model = "rateModel"), function(x,model,scale=1,threads=1) {
  return(model@phylogeny$ll(model@alleleData$alleleData@data@x,
                              model@rateDM@x,
                              model@piDM@x,scale,threads))
})

#' @name scaledLL
#' @rdname scaledLL
methods::setMethod("scaledLL", signature(x="numeric",model = "rateModel"), function(x,model,scale=1,threads=1) {
  setParams(model,x,which(!model@fixed)-1)
  return(model@phylogeny$ll(model@alleleData$alleleData@data@x,
                          model@rateDM@x,
                          model@piDM@x,scale,threads))
})

#' Compute log-likelihood
#'
#' Compute the log-likelihood for given rate model object.
#' @param x a set of parameters to compute the log-likelihood for (optional)
#' @param model rateModel object
#' @param threads number of threads to use
#' @return scaled log-likelihood
#' @name ll
#' @include rateModel-class.R
#' @rdname ll
#' @export
methods::setGeneric("ll", function(x,model,threads=1,...) {
  standardGeneric("ll")
})


#' @name ll
#' @rdname ll
methods::setMethod("ll", signature(x="missing",model = "rateModel"), function(x,model,threads=1) {
  return(model@phylogeny$ll(model@alleleData$alleleData@data@x,
                            model@rateDM@x,
                            model@piDM@x,1,threads))
})

#' @name ll
#' @rdname ll
methods::setMethod("ll", signature(x="numeric",model = "rateModel"), function(x,model,threads=1) {
  setParams(model,x,which(!model@fixed)-1)
  return(model@phylogeny$ll(model@alleleData$alleleData@data@x,
                            model@rateDM@x,
                            model@piDM@x,1,threads))
})
