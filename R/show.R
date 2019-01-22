#' Replace default print format
#'
#' Method to replace default print appearance for 
#' @param x alleleData or rateModel object
#' @name show
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @return A tree of class \code{phylo}
#' @rdname show

#' @name show
#' @rdname show
methods::setMethod("show", "alleleData", function(object) {
  print(str(object))
})

#' @name show
#' @rdname show
methods::setMethod("show",  "rateModel", function(object) {
  print(str(object))
})