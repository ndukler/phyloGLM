#' Node marginals
#'
#' Compute the posterior marginal distribution across alleles per node
#' @param model rateModel object
#' @param threads number of threads to use
#' @return a list of per node marginal posteriors
#' @name marginal
#' @include rateModel-class.R
#' @rdname marginal
#' @export
methods::setGeneric("marginal", function(model,threads=1,...) {
  standardGeneric("marginal")
})

#' @name marginal
#' @rdname marginal
methods::setMethod("marginal", signature(model = "rateModel"), function(model,threads=1) {
  return(model@phylogeny$marginal(model@alleleData$alleleData@data@x,
                              model@rateDM@x,
                              model@piDM@x,threads))
})