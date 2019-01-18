methods::setGeneric("getAlleleData", function(obj) {
  standardGeneric("getAlleleData")
})

#' Returns alleleData
#'
#' Gets the alleleData from a rate model object
#' @param obj rateModel 
#' @name getAlleleData
#' @return alleleData object
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("getAlleleData", signature(obj = "rateModel"), function(obj) {
  return(obj@alleleData$alleleData)
})

