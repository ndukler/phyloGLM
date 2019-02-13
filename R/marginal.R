#' Node marginals
#'
#' Compute the posterior marginal distribution across alleles per node
#' @param obj rateModel object
#' @param threads number of threads to use
#' @return a list of per node marginal posteriors
#' @name marginal
#' @include rateModel-class.R
#' @rdname marginal
#' @export
methods::setGeneric("marginal", function(obj,threads=1,...) {
  standardGeneric("marginal")
})

#' @name marginal
#' @rdname marginal
methods::setMethod("marginal", signature(obj = "rateModel"), function(obj,threads=1) {
  return(obj@phylogeny$marginal(obj@alleleData$alleleData@data@x,
                              obj@rateDM@x,
                              obj@piDM@x,threads))
})