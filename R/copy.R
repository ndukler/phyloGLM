#' Copy rateModel
#'
#' Creates a copy of a rate model
#' @param model a rateModel object
#' @name copy
#' @include rateModel-class.R
#' @rdname copy
#' @return a copyed of a rateModel object
#' @examples
#' 
#' @export
methods::setGeneric("copy", function(model) {
  standardGeneric("copy")
})

#' @name copy
#' @rdname copy
methods::setMethod("copy", signature(model = "rateModel"), function(model) {
  return(unpack(pack(model)))
})
