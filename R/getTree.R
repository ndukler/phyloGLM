#' Gets tree
#'
#' Returns phylo object
#' @param model alleleData or rateModel object
#' @name getTree
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @return A tree of class \code{phylo}
#' @rdname getTree
#' @examples
#' 
#' @export
methods::setGeneric("getTree", function(model) {
  standardGeneric("getTree")
})

#' @name getTree
#' @rdname getTree
methods::setMethod("getTree", signature(model = "alleleData"), function(model) {
  return(model@tree)
})

#' @name getTree
#' @rdname getTree
methods::setMethod("getTree", signature(model = "rateModel"), function(model) {
  return(getTree(getAlleleData(model)))
})
