#' Extract site information table
#'
#' Extracts site information table
#' @docType methods
#' @param obj alleleData or rateModel object
#' @rdname getSiteInfo
#' @return data.table containing site info
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @export
methods::setGeneric("getSiteInfo", function(obj) {
  standardGeneric("getSiteInfo")
})

#' @rdname getSiteInfo
#' @aliases getSiteInfo,alleleData,alleleData-method
methods::setMethod("getSiteInfo", signature(obj = "alleleData"), function(obj) {
  return(obj@siteInfo)
})

#' @rdname getSiteInfo
#' @aliases getSiteInfo,rateModel,rateModel-method
methods::setMethod("getSiteInfo", signature(obj = "rateModel"), function(obj) {
  return(getSiteInfo(getAlleleData(obj)))
})