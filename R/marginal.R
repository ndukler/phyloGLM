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
methods::setGeneric("marginal", function(model, threads = 1, ...) {
  standardGeneric("marginal")
})

#' @name marginal
#' @rdname marginal
methods::setMethod("marginal", signature(model = "rateModel"), function(model, threads = 1) {
  return(model@phylogeny$marginal(
    model@alleleData$alleleData@data@x,
    model@rateDM@x,
    model@piDM@x, threads
  ))
})

#' Expected transitions per site
#'
#' Compute the expected number of transitions aggregated across either nodes or edges
#' @param model rateModel object
#' @param aggregate How expected numbers of transitions should be aggregated (either by edge or by node)
#' @param threads number of threads to use
#' @return a list of per node/edge matricies containing the expected number of allelic transitions. The
#' rows are the ancestral allele and the columns are the descendant allele
#' @name marginalTransitions
#' @include rateModel-class.R
#' @rdname marginalTransitions
#' @export
methods::setGeneric("marginalTransitions", function(model, aggregate = c("edge", "node"), threads = 1) {
  standardGeneric("marginalTransitions")
})

#' @name marginalTransitions
#' @rdname marginalTransitions
methods::setMethod("marginalTransitions", signature(model = "rateModel"), function(model, aggregate = c("edge", "node"),
                                                                                   threads = 1) {
  ## Set default for aggregations
  aggregate <- aggregate[1]
  if (aggregate == "edge") {
    et <- model@phylogeny$edgewiseMarginalTransitions(
      model@alleleData$alleleData@data@x,
      model@rateDM@x,
      model@piDM@x, threads
    )
  } else if (aggregate == "node") {
    et <- model@phylogeny$nodewiseMarginalTransitions(
      model@alleleData$alleleData@data@x,
      model@rateDM@x,
      model@piDM@x, threads
    )
  } else {
    stop("Invalid aggregation option")
  }
  ## Collapse transition counts to list of matricies
  out <- lapply(et, function(x) do.call(what = "rbind", args = x))
  return(out)
})
