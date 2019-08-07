## Send request to ape package people to do this rather than doing it myself
## setOldClass("phylo")

#' Class alleleData
#'
#' Class \code{alleleData} holds the allelic probabilities and corresponding phylogenetic tree
#'
#' @name alleleData-class
#' @rdname alleleData-class
#' @include alleleDataValidityCheck.R
#' @include stlMatrix-class.R
#' @importClassesFrom data.table data.table
#' @exportClass alleleData
methods::setClass("alleleData",
  slots = c(data = "stlMatrix", tree = "ANY", siteInfo = "data.table", nAlleles = "numeric", nSpecies = "numeric", nSites = "numeric"),
  validity = alleleDataValidityCheck
)

#' alleleData
#'
#' Contructs an object that holds the allele data and the phylogenetic tree
#' @param data A list of matricies, where each element in the list is a species, each column in a matrix is an allele, and each row is a site
#' @param tree The tree corresponding the the list of allele data passed in. Tip names must match list names.
#' @param logProb TRUE/FALSE indicating whether the passed probabilities are in log space (default=FALSE)
#' @param siteInfo A data.frame containing classification information about each site, where each column characteristic, and each site is a row. Defaults to assuming all sites have the same info.
#' @name alleleData
#' @return an alleleData object
#' @examples
#' 
#' @export
alleleData <- function(data, tree, siteInfo = NULL, logProb = FALSE) {
  ## **Some checks prior to beginning construction**
  ## Check that tree is phylo type
  if (class(tree) != "phylo") {
    stop("Tree must be of class phylo")
  }
  ## Check that there are the same number of tree tips as there are list elements
  if (length(data) != sum(!tree$edge[, 2] %in% tree$edge[, 1])) {
    stop("Different numbers of tree tips and list elements in data")
  }
  ## Check that the labels at the tips of the tree match the list names
  if (!all(names(data) %in% tree$tip.label) && !is.null(names(data)) && !is.null(tree$tip.label)) {
    stop("tree$tip.label and names(data) do not contain the same species")
  }
  ## Check that there are the same number of rows in all data entries
  if (length(unique(unlist(lapply(data, nrow)))) != 1) {
    stop("Differing numbers of sites for different species in data")
  }
  ## Check that siteInfo is a data.frame (or data.table)
  if (!is.null(siteInfo) && !is.data.frame(siteInfo)) {
    stop("siteInfo must be a data.frame or a data.table")
  }

  
  ## Post-order tree
  tree <- ape::reorder.phylo(tree, "postorder")

  ## Compute number of alleles and species
  nAlleles <- ncol(data[[1]])
  nSpecies <- length(data)
  ## Convert list to matrix
  dataMatrix <- do.call("cbind", data[tree$tip.label])
  colnames(dataMatrix) <- paste0(
    rep(tree$tip.label, unlist(lapply(data[tree$tip.label], ncol))), ".",
    unlist(lapply(data[tree$tip.label], function(x) 1:ncol(x)))
  )
  ## If site info is NULL, create faux site info that groups all sites together
  if (is.null(siteInfo)) {
    siteInfo <- data.table::data.table(default = rep(1, nrow(dataMatrix)))
  }
  
  ## Ensure that all character columns have been converted to factors
  siteInfo <- data.table::data.table(siteInfo, stringsAsFactors = T)
  
  ## Log data if necessary
  if (!logProb) {
    dataMatrix <- log(dataMatrix)
  }

  ##
  ## Do data validity checks
  ##

  ## Check that all numbers are valid probabilities
  if (any(dataMatrix > 0) || any(is.na(dataMatrix)) || any(is.nan(dataMatrix))) {
    stop("Non-valid log-proababilities in alleleData@data. All values must be in [-Inf,0].")
  }
  ## Check that there are the same number of alleles in each species
  if (table(unique(table(gsub("\\.\\d+$", "", x = colnames(dataMatrix))))) != 1) {
    stop("Differing numbers of alleles between species")
  }
  ## Check that the labels at the tips of the tree match the transformed names in the matrix
  if (any(unique(gsub("\\.\\d+$", "", x = colnames(dataMatrix))) != tree$tip.label)) {
    stop("Species names do not match in the tree and the data matrix")
  }
  ## Check that there are the expected number of columns in the data matrix
  if (ncol(dataMatrix) != (nAlleles * nSpecies)) {
    stop("The number of columns in the data matrix is not the expected nSpecies*nAlleles")
  }
  
  ## Check that there are the same number of rows in the covarirate and data matricies
  if(nrow(dataMatrix)!=nrow(siteInfo)){
    stop("Differing numbers of sites (rows) in data and siteInfo")
  }
  
  ## End checks

  ## Get nSites
  nSites <- nrow(dataMatrix)

  ## Build object
  methods::new("alleleData",
    data = stlMatrix(dataMatrix), tree = tree, siteInfo = data.table::as.data.table(siteInfo),
    nSpecies = nSpecies, nAlleles = nAlleles, nSites = nSites
  )
}
