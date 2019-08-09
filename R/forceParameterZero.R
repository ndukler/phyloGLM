#' Remove parameter from model
#'
#' Sets a parameter to 0 and sets its value to be fixed so it will not be fit
#' @param obj rateModel object
#' @param index a 0-based integer vector of parameter indicies to be changed (note: parameters are 0 indexed)
#' @param verbose if TRUE prints out information about variable being frozen
#' @rdname forceParameterZero
#' @name forceParameterZero
#' @include rateModel-class.R
#' @export
methods::setGeneric("forceParameterZero", function(obj,index, verbose = TRUE) {
  standardGeneric("forceParameterZero")
})

#' @name setParams
#' @rdname setParams
#' @aliases setParams,rateModel,rateModel-method
methods::setMethod("forceParameterZero", signature(obj = "rateModel"), function(obj, index, verbose = TRUE) {
  updateFixed(old = obj@fixed, replacement = TRUE, index = index)
  setParams(obj = obj, x = 0, index = index)
  if(verbose){
    rateIndex = obj@phylogeny$getRateIndex()
    piIndex = obj@phylogeny$getPiIndex()
    if(index %in% rateIndex$idx){
      info = c(type = "rate",rateIndex[rateIndex$idx==index,])  
    } else {
      info = c(type = "pi",rateIndex[piIndex$idx==index,])  
    }
  }
  if(verbose)
    with(info,print(paste0("Variable <",name,"> of type <",type,"> in group <",group,"> has been fixed to 0.")))
  
})
