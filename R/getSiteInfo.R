#' Extract site information table
#'
#' Extracts site information table
#' @param model alleleData or rateModel object
#' @name getSiteInfo
#' @rdname getSiteInfo-methods
#' @return data.table containing site info
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @export
methods::setGeneric("getSiteInfo", function(model) {
  standardGeneric("getSiteInfo")
})

#' @name getSiteInfo
#' @rdname getSiteInfo-methods
#' @aliases getSiteInfo,alleleData,alleleData-method
methods::setMethod("getSiteInfo", signature(model = "alleleData"), function(model) {
  return(model@siteInfo)
})

#' @name getSiteInfo
#' @rdname getSiteInfo-methods
#' @aliases getSiteInfo,rateModel,rateModel-method
methods::setMethod("getSiteInfo", signature(model = "rateModel"), function(model) {
  return(getSiteInfo(getAlleleData(model)))
})
