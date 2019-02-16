#' Compute sitewise log-likelihoods
#'
#' Compute sitewise log-likelihoods
#' @param model rateModel object
#' @param threads number of threads to use
#' @rdname siteLL
#' @name siteLL
#' @return vector of parameter values
#' @include rateModel-class.R
#' @export 
methods::setGeneric("siteLL", function(model,threads=1) {
  standardGeneric("siteLL")
})

#' @name siteLL
#' @rdname siteLL
#' @aliases siteLL,rateModel,rateModel-method 
methods::setMethod("siteLL", signature(model = "rateModel"), function(model,threads=1) {
  # if(is.null(sites)){
  #   sites=1:nrow(model@rateDM)
  # } else if (any(sites==0)){
  #   stop("sites is 1-indexed (at least for now)")
  # }
  return(model@phylogeny$siteLL(model@alleleData$alleleData@data@x,
                              model@rateDM@x,
                              model@piDM@x,threads))
})
