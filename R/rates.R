#' Compute rates
#'
#' Compute rates of current model across specified sites
#' @param model rateModel object
#' @param sites integer vector of site indicies
#' @return table of rates along with the corresponding edgeGroup and site
#' @name rates
#' @include rateModel-class.R
#' @rdname rates
#' @export
methods::setGeneric("rates", function(model,sites=NULL) {
  standardGeneric("rates")
})


#' @name rates
#' @rdname rates
methods::setMethod("rates", signature(model = "rateModel"), function(model,sites=NULL) {
  ## If no sites are selected, default to all of them
  if(is.null(sites)){
    sites=1:model@alleleData$alleleData@data@nrow
  }
  if(any(sites<1 | sites > model@alleleData$alleleData@data@nrow)){
    stop("Out of range site specified")
  }
  ## Select a representitive child for each edge group and adjust index
  repChildren=model@edgeGroups[!duplicated(model@edgeGroups$edgeGroup),.(child,edgeGroup)]
  ##compute rates
  out=data.table::as.data.table(expand.grid(edgeGroup=as.integer(repChildren$edgeGroup),site=as.integer(sites)))
  out=merge(repChildren,out,by="edgeGroup")
  out[,rate:=model@phylogeny$rate(child-1,model@rateDM[site,]),by=c("child","sites")]
  out[,child:=NULL]
  return(out)
})