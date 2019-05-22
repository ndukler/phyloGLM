#' Compute site-wise stationary distribution
#'
#' Compute the site-wise allelic stationary distribution across specified sites
#' @param model rateModel object
#' @param sites integer vector of site indicies
#' @return table of allelic probabilities along with the corresponding allele and site
#' @name pi
#' @include rateModel-class.R
#' @rdname pi
#' @export
methods::setGeneric("pi", function(model, sites = NULL) {
  standardGeneric("pi")
})


#' @name pi
#' @rdname pi
methods::setMethod("pi", signature(model = "rateModel"), function(model, sites = NULL) {
  ## If no sites are selected, default to all of them
  if (is.null(sites)) {
    sites <- 1:model@alleleData$alleleData@data@nrow
  }
  if (any(sites < 1 | sites > model@alleleData$alleleData@data@nrow)) {
    stop("Out of range site specified")
  }
  ## compute pi
  ## Compute stationary distribution
  out_pi <- melt(data.table::data.table(model@phylogeny$piV_Rcpp(sites - 1, model@piDM@x))[,site:=sites],id.vars="site")
  data.table::setnames(out_pi,c("variable","value"),c("allele","probability"))
  out_pi[,allele:=as.integer(gsub("V","",allele))]
  setkeyv(out_pi,cols = c("site","allele"))
  return(out)
})
