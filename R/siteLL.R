#' Compute sitewise log-likelihoods
#'
#' Compute sitewise log-likelihoods
#' @param x rateModel object
#' @param sites integer vectore of sites that log-likelihoods should be computed for (note: 1-indexed)
#' @rdname siteLL
#' @name siteLL
#' @return vector of parameter values
#' @include rateModel-class.R
#' @export 
methods::setGeneric("siteLL", function(obj,sites=NULL,threads=1) {
  standardGeneric("siteLL")
})

#' @name siteLL
#' @rdname siteLL
#' @aliases siteLL,rateModel,rateModel-method 
methods::setMethod("siteLL", signature(obj = "rateModel"), function(obj,sites=NULL,threads=1) {
  if(is.null(sites)){
    sites=1:nrow(obj@rateDM)
  } else if (any(sites==0)){
    stop("sites is 1-indexed (at least for now)")
  }
  return(obj@phylogeny$siteLL(obj@alleleData$alleleData@data[sites,,drop=FALSE],
                              obj@rateDM[sites,,drop=FALSE],
                              obj@piDM[sites,,drop=FALSE],threads))
})
