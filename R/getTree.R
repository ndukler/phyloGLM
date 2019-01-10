#' Gets tree
#'
#' Returns phylo object
#' @param obj alleleData or rateModel object
#' @name getTree
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @return A tree of class \code{phylo}
#' @rdname getTree
#' @examples
#' 
#' @export
methods::setGeneric("getTree", function(obj) {
  standardGeneric("getTree")
})

#' @name getTree
#' @rdname getTree
methods::setMethod("getTree", signature(obj = "alleleData"), function(obj) {
  return(obj@tree)
})

#' @name getTree
#' @rdname getTree
methods::setMethod("getTree", signature(obj = "rateModel"), function(obj) {
  return(getTree(getAlleleData(obj)))
})