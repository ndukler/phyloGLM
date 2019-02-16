#' Get edge table
#'
#' Returns edge table.
#' @param model alleleData or rateModel object
#' @rdname getEdgeTable
#' @name getEdgeTable
#' @return data.table containing edge ids
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @export 
methods::setGeneric("getEdgeTable", function(model) {
  standardGeneric("getEdgeTable")
})

#' @name getEdgeTable
#' @rdname getEdgeTable
#' @aliases getEdgeTable,alleleData,alleleData-method
methods::setMethod("getEdgeTable", signature(model = "alleleData"), function(model) {
  data.table::data.table(parent=model@tree$edge[,1],child=model@tree$edge[,2])
})

#' @name getEdgeTable
#' @rdname getEdgeTable
#' @aliases getEdgeTable,rateModel,rateModel-method 
methods::setMethod("getEdgeTable", signature(model = "rateModel"), function(model) {
  return(model@edgeGroups)
})