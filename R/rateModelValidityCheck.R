#' @description  Validity checker for alleleData-class
rateModelValidityCheck <- function(object) {
  errors <- c()
  if (!setequal(names(object@alleleData), c("alleleData"))) {
    errors <- c(errors, "object@alleleData must have one binding, \'alleleData\'")
  } else if (class(object@alleleData$alleleData) != "alleleData") {
    errors <- c(errors, "object@alleleData$alleleData is not of class allele Data")
  }
  if (!setequal(colnames(object@edgeGroups), c("parent", "child", "edgeGroup"))) {
    errors <- c(errors, "Colnames of edgeGroups must be \'parent\', \'child\',\'edgeGroup\'")
  }
  if (!setequal(data.table::key(object@edgeGroups), c("child"))) {
    errors <- c(errors, "Key of edgeGroups must be \'child\'")
  }
  if (class(object@phylogeny) != "Rcpp_phylogeny") {
    errors <- c(errors, "phylogeny must be a \`Rcpp_phylogeny\` object")
  }
  ## Lock the environment and all the bindings in it if all tests passed
  if (length(errors) == 0) {
    lockEnvironment(env = object@alleleData, bindings = TRUE)
  }
  ## Return errors if there are any
  if (length(errors) == 0) TRUE else errors
}
