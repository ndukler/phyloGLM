#' @description  Validity checker for alleleData-class
alleleDataValidityCheck <- function(object){
  errors=c()
  if(class(object@tree) != "phylo"){
    errors= c(errors,"Tree must be of class phylo")
  }
  ## Check that there are the same number of alleles in each species
  if(table(unique(table(gsub("\\.\\d+$", "",x=colnames(object@data)))))!=1){
    errors=c(errors,"Differing numbers of alleles between species")
  }
  ## Check that there are a minimum of two alleles
  if(object@nAlleles <= 1){
    errors=c(errors,"Less than two alleles per species")
  }
  ## Check that there are the expected number of columns in the data matrix
  if(ncol(object@data)!=(object@nAlleles*object@nSpecies)){
    errors=c(errors,"The number of columns in the data matrix is not the expected nSpecies*nAlleles") 
  }
  ## Check that the labels at the tips of the tree match the transformed names in the matrix
  if(any(unique(gsub("\\.\\d+$", "",x=colnames(object@data))) != object@tree$tip.label)){
    errors=c(errors,"Species names do not match in the tree and the data matrix") 
  }
  ## Check that all numbers are valid probabilities
  if(any(object@data > 0) || any(is.na(object@data)) || any(is.nan(object@data))){
    errors=c(errors,c("Non-valid log-proababilities in alleleData@data. All values must be in [-Inf,0]."))
  }
  ## Check that tree is unrooted
  # if(ape::is.rooted(object@tree)){
  #   errors=c(errors,c("Tree must be unrooted."))  
  # }
  ## Check that tree is post-ordered
  if(any(object@tree$edge != ape::reorder.phylo(object@tree, "postorder")$edge)){
    errors=c(errors,c("Tree must be sorted in post-order travel order."))  
  }
  ## Check that the number of sites are correct
  if(nrow(object@data)!=object@nSites){
    errors=c(errors,c("Incorrect number of sites"))  
  }
  ## Check that there are the same number of rows in the siteInfo data.frame as there are rows in the data
  if (length(errors) == 0) TRUE else errors
}