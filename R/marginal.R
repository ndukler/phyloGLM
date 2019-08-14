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
#' @param aggregate How expected numbers of transitions should be aggregated (within edge or node)
#' @param excludeNodes Indicates nodes that should be excluded from the per-edge aggregated statistic 
#' (1-based). Only usable with aggregate = "edge".
#' @param excludeEdges Indicates edges that should be excluded from the per-node aggregated statistic 
#' Only usable with aggregate = "node". (edge ids based on integer id of child, use plotTree(model) 
#' to view).
#' @param threads number of threads to use
#' @return a list of per node/edge matricies containing the expected number of allelic transitions. The
#' rows are the ancestral allele and the columns are the descendant allele. The list ids are a character
#' string that are either a site or edge (child) id depending on how you aggregated
#' @name marginalTransitions
#' @include rateModel-class.R
#' @rdname marginalTransitions
#' @export
methods::setGeneric("marginalTransitions", function(model, aggregate = c("edge", "node"), excludeNodes = integer(0),
                                                    excludeEdges = integer(0), threads = 1) {
  standardGeneric("marginalTransitions")
})

#' @name marginalTransitions
#' @rdname marginalTransitions
methods::setMethod("marginalTransitions", signature(model = "rateModel"), function(model, aggregate = c("edge", "node"),
                                                                                   excludeNodes,excludeEdges, 
                                                                                   threads = 1) {
  ## Set default for aggregations
  aggregate <- aggregate[1]
  if (aggregate == "edge") {
    ## Check that all elements of exclude are in range
    if(any(excludeNodes > model@rateDM@nrow) | any(excludeNodes < 1)){
      stop("Invalid sites specified in exclude (less than 1 or greater than num. sites).")  
    }
    if(length(excludeEdges)>0){
      warning("excludeEdges not implemented when aggregate = \'edge\' and will be ignored.")
    }
    et <- model@phylogeny$edgewiseMarginalTransitions(
      model@alleleData$alleleData@data@x,
      model@rateDM@x,
      model@piDM@x,
      excludeNodes-as.integer(1),
      threads
    )
    
    ## Collapse transition counts to list of matricies
    out <- lapply(et, function(x) do.call(what = "rbind", args = x))
    
    ## Add as characters the ids of the edges as list names
    labels <- as.character(getTree(model)$edge[,2])
  } else if (aggregate == "node") {
    ## Check that all elements of exclude are in range
    if(any(!excludeEdges %in% getTree(model)$edge[,2])){
      stop("Invalid edged specified in exclude. To see valid ids run getTree(rate_model)$edge[,2] or 
           plotTree(rate_model).")  
    }
    if(length(excludeNodes)>0){
      warning("excludeNodes not implemented when aggregate = \'node\' and will be ignored.")
    }
    et <- model@phylogeny$nodewiseMarginalTransitions(
      model@alleleData$alleleData@data@x,
      model@rateDM@x,
      model@piDM@x, 
      excludeEdges-as.integer(1),
      threads
    )
    
    ## Collapse transition counts to list of matricies
    out <- lapply(et, function(x) do.call(what = "rbind", args = x))
    
    ## Add as characters the ids of the nodes as list names
    labels <- as.character(1:model@rateDM@nrow)

  } else {
    stop("Invalid aggregation option")
  }
  names(out) <- labels
  return(out)
})
