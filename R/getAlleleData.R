methods::setGeneric("getAlleleData", function(model) {
  standardGeneric("getAlleleData")
})

#' Returns alleleData
#'
#' Gets the alleleData from a rate model object
#' @param model rateModel 
#' @name getAlleleData
#' @return alleleData object
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("getAlleleData", signature(model = "rateModel"), function(model) {
  return(model@alleleData$alleleData)
})

