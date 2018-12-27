#' Get edge table
#'
#' Returns edge table.
#' @docType methods
#' @param x alleleData or rateModel object
#' @rdname getEdgeTable
#' @return data.table containing edge ids
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @export
methods::setGeneric("getEdgeTable", function(obj) {
  standardGeneric("getEdgeTable")
})

#' @rdname getEdgeTable
#' @aliases getEdgeTable,alleleData,alleleData-method
methods::setMethod("getEdgeTable", signature(obj = "alleleData"), function(obj) {
  data.table::data.table(parent=obj@tree$edge[,1],child=obj@tree$edge[,2])
})

#' @rdname getEdgeTable
#' @aliases getEdgeTable,rateModel,rateModel-method 
methods::setMethod("getEdgeTable", signature(obj = "rateModel"), function(obj) {
  return(obj@edgeGroups)
})