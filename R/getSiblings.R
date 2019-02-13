#' Get siblings
#'
#' Get the siblings of a given node on a tree
#' @param tr phylo object
#' @param node the node number
#' @return The node numbers of sibbling nodes
#' @name getSiblings
#' @include rateModel-class.R
#' @rdname getSiblings
#' @export
getSiblings = function(tr,node){
  if(!inherits(tr,"phylo")) {
    stop("tree should be an object of class \"phylo\".")
  }
  parent=tr$edge[tr$edge[,2]==node,1]
  return(setdiff(tr$edge[tr$edge[,1]==parent,2],node))
}