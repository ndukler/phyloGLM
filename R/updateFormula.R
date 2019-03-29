#' Update formula
#'
#' Updates the rate and/or pi formula of a rate model
#' @param model a rateModel object
#' @param rateFormula a new rate formula
#' @param piFormula a new pi formula
#' @name updateFormula
#' @include rateModel-class.R
#' @rdname updateFormula
#' @return a copy of the rateModel object with updated formulas (matching parameters carry over values)
#' @examples
#' 
#' @export
methods::setGeneric("updateFormula", function(model, rateFormula, piFormula) {
  standardGeneric("updateFormula")
})

#' @name updateFormula
#' @rdname updateFormula
methods::setMethod("updateFormula", signature(model = "rateModel"), function(model, rateFormula = NULL, piFormula = NULL) {
  ## Serialize model for easy copying
  pm <- pack(model)
  ## Extract parameter indicies
  oldRateIndex <- data.table::as.data.table(model@phylogeny$getRateIndex())[, value := pm@params[idx + 1]]
  oldPiIndex <- data.table::as.data.table(model@phylogeny$getPiIndex())[, value := pm@params[idx + 1]]
  ## Update formulas or use existing ones as defaults
  if (is.null(rateFormula)) {
    rateFormula <- pm@rateFormula
  }
  if (is.null(piFormula)) {
    piFormula <- pm@piFormula
  }
  ## Create new rate model object with updated formulas
  r <- rateModel(
    data = alleleData(pm@data, pm@tree, pm@siteInfo, logProb = TRUE), rateFormula = rateFormula,
    piFormula = piFormula, lineageTable = pm@lineageTable, rateBounds = pm@rateBounds
  )
  ## Get new parameter indicies
  newRateIndex <- data.table::as.data.table(r@phylogeny$getRateIndex())
  newPiIndex <- data.table::as.data.table(r@phylogeny$getPiIndex())
  ## Merge new and old parameter indicies on name and group
  mergedRateIndex <- merge(newRateIndex, oldRateIndex[, .(name, group, value)], by = c("name", "group"))
  mergedPiIndex <- merge(newPiIndex, oldPiIndex[, .(name, group, value)], by = c("name", "group"))
  ## Set parameters to matched values
  setParams(obj = r, x = mergedRateIndex$value, index = mergedRateIndex$idx)
  setParams(obj = r, x = mergedPiIndex$value, index = mergedPiIndex$idx)
  return(r)
})
