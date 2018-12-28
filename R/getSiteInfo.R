#' Extract site information table
#'
#' Extracts site information table
#' @param obj alleleData or rateModel object
#' @name getSiteInfo
#' @rdname getSiteInfo-methods
#' @return data.table containing site info
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @export
methods::setGeneric("getSiteInfo", function(obj) {
  standardGeneric("getSiteInfo")
})

#' @name getSiteInfo
#' @rdname getSiteInfo-methods
#' @aliases getSiteInfo,alleleData,alleleData-method
methods::setMethod("getSiteInfo", signature(obj = "alleleData"), function(obj) {
  return(obj@siteInfo)
})

#' @name getSiteInfo
#' @rdname getSiteInfo-methods
#' @aliases getSiteInfo,rateModel,rateModel-method
methods::setMethod("getSiteInfo", signature(obj = "rateModel"), function(obj) {
  return(getSiteInfo(getAlleleData(obj)))
})