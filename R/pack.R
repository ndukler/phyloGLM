#' Class packedModel
#'
#' Class \code{packedModel} holds all the data necessary to construct a rate model
#'
#' @slot alleleData an allele data object in a locked environment
#' @slot edgeGroups a data.table with four columns: parent, child, edgeGroup. All parent-child combos must be valid for alleleData object.
#' @slot rateFormula A formula that uses the variables in the allele data object to compute turnover rate
#' @slot piFormula A formula that uses the variables in the allele data object to compute pi
#' @slot rateDM the design matrix for the rate
#' @slot piDM the design matrix for pi
#' @slot phylogeny parameter object holding all the parameters for phylogenetic computations
#' @slot fixed a logical vector indicating which variables are fixed
#'
#' @name packedModel-class
#' @rdname packedModel-class
#' @importClassesFrom data.table data.table
#' @exportClass packedModel
methods::setClass("packedModel", slots = c(
  tree = "ANY", data = "list",
  siteInfo = "data.table", rateFormula = "formula", piFormula = "formula",
  lineageTable = "data.table", params = "numeric", fixed = "logical",
  rateBounds = "numeric"
))

#' Pack rateModel
#'
#' Packs a rateModel into an S4 object that can be saved
#' @param model a rateModel object
#' @name pack
#' @include rateModel-class.R
#' @rdname pack
#' @return a packedModel object
#' @examples
#' 
#' @export
methods::setGeneric("pack", function(model) {
  standardGeneric("pack")
})

#' @name pack
#' @rdname pack
methods::setMethod("pack", signature(model = "rateModel"), function(model) {
  ## Extract all the objects needed to rebuild the alleleData object
  ad <- getAlleleData(model)
  tree <- getTree(ad)
  dat <- ad@data[1:ad@data@nrow,,drop=FALSE]
  a <- rep(tree$tip.label, each = ad@nAlleles)
  data <- lapply(split(seq_along(a), a)[tree$tip.label], function(ind, m) m[, ind, drop=FALSE], m = dat)
  rm(dat)
  siteInfo <- ad@siteInfo
  ## Extract all the objects needed to rebuild the rateModel object
  rateFormula <- model@rateFormula
  piFormula <- model@piFormula
  lineageTable <- model@edgeGroups
  params <- getParams(model)
  fixed <- model@fixed
  rateBounds <- getRateBounds(model)
  return(new("packedModel",
    tree = tree, data = data, siteInfo = siteInfo, rateFormula = rateFormula, piFormula = piFormula,
    lineageTable = lineageTable, params = params, fixed = fixed, rateBounds = rateBounds
  ))
})

#' Unpack rateModel
#'
#' Unpacks a packedModel in to a rateModel
#' @param pm a packedModel object
#' @name unpack
#' @include rateModel-class.R
#' @rdname unpack
#' @return a rateModel object
#' @examples
#' 
#' @export
methods::setGeneric("unpack", function(pm) {
  standardGeneric("unpack")
})

#' @name unpack
#' @rdname unpack
methods::setMethod("unpack", signature(pm = "packedModel"), function(pm) {
  ## Build a rate model object
  r <- rateModel(
    data = alleleData(pm@data, pm@tree, pm@siteInfo, logProb = TRUE), rateFormula = pm@rateFormula,
    piFormula = pm@piFormula, lineageTable = pm@lineageTable, rateBounds = pm@rateBounds
  )
  setParams(r, pm@params, (1:length(pm@params)) - 1)
  setFixed(r,pm@fixed, 0:(length(pm@fixed)-1))
  return(r)
})
