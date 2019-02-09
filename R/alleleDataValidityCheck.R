#' @description  Validity checker for alleleData-class
alleleDataValidityCheck <- function(object){
  errors=c()
  if(class(object@tree) != "phylo"){
    errors= c(errors,"Tree must be of class phylo")
  }
  ## Check that there are a minimum of two alleles
  if(object@nAlleles <= 1){
    errors=c(errors,"Less than two alleles per species")
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
  # if(nrow(object@data)!=object@nSites){
  #   errors=c(errors,c("Incorrect number of sites"))  
  # }
  ## Check that there are the same number of rows in the siteInfo data.frame as there are rows in the data
  if (length(errors) == 0) TRUE else errors
}