#' Simulate data
#'
#' Simulates data without an error model
#' @param tr tree
#' @param covariateTable a data.frame where each column is a covariate and each row 
#' is a site (columns have covariate names)
#' @param rateFormula a formula object specifying how covariates are combined to compute the rate
#' @param rateParams a matrix of parameters where the columns correspond to the covariate label and the rows are
#' the edgeGroup from the lineage table
#' @param piFormula a formula object specifying how covariates are combined to compute the stationary distribution
#' @param piParams a matrix of parameters where the columns correspond to the covariate label and the rows are
#' each seperate alleles (note nAlleles = nrow(piParams)+1)
#' @param lineageTable a data.frame with three columns, parent, child, and edgeGroup. Can be used to specify different
#' rate parameters for different edges
#' @name simulateSites
#' @return list with the simulated data and the number of alleles
#' @export

rateFormula=formula(~x1*x2+0)
covariateTable=data.frame(x1=sample(c(0,1),size = 1000,replace = T),x2=sample(c(0,1),size = 1000,replace = T))

simulateSites <- function(tr,covariateTable,rateFormula,rateParams=NULL,piFormula=NULL,piParams=NULL,lineageTable=NULL){
  ## **Start parameter tests**
  
  ## Check if piFormula is specified. If not, set to same as rateFormula
  if(is.null(piFormula)){
    piFormula=rateFormula
  }
  
  ## Create default lineageTable and perform checks if it is already specified
  if(is.null(lineageTable)){
    lineageTable=data.table::data.table(parent=tr$edge[,1],child=tr$edge[,2],edgeGroup=1)
  } else if (!is.data.frame(lineageTable) || colnames(lineageTable) != c("parent","child","edgeGroup")) {
    stop("lineageTable must be a data.frame with colnames: parent, child, edgeGroup")
  } else if (!is.integer(lineageTable$edgeGroup)){
    stop("edgeGroup must be integer valued.")
  } else if (nrow(lineageTable) != nrow(tr$edge)){
    stop("lineageTable and the tree must contain the same number of edges")
  } else if (any(as.matrix(lineageTable[order(child),.(parent,child)])!=tr$edge[order(tr$edge[,2]),])){
    stop("lineageTable and the specified tree have one or more different parent-child edges.")
  }
  
  ## Warn user if the formula contains an intercept and provide info on how to avoid 
  if(attr(terms(rateFormula),"intercept")!=0){
    warning("rateFormula has an intercept term. To prevent this specify formula with +0 at the end.")
  }
  if(attr(terms(piFormula),"intercept")!=0){
    warning("piFormula has an intercept term. To prevent this specify formula with +0 at the end.")
  }
  ## Check that the covariate table has all variables required by the formulas
  if(!all(all.vars(rateFormula) %in% colnames(covariateTable))){
    missing=setdiff(all.vars(rateFormula),colnames(covariateTable))
    stop(paste("Covariates in the rate formula missing from the covariate table:",paste(missing,collapse = ",")))
  }
  if(!all(all.vars(piFormula) %in% colnames(covariateTable))){
    missing=setdiff(all.vars(piFormula),colnames(covariateTable))
    stop(paste("Covariates in the pi formula missing from the covariate table:",paste(missing,collapse = ",")))
  }
  
  ## Check that the form of the parameter matricies are correct
  nFeaturesRate=length(attr(terms(rateFormula),"term.labels")) + attr(terms(rateFormula),"intercept")
  nFeaturesPi=length(attr(terms(piFormula),"term.labels")) + attr(terms(piFormula),"intercept")
  nEdgeGroup=length(unique(lineageTable$edgeGroup))
  if(is.null(rateParams)){
    rateParams=matrix(1,ncol = nFeaturesPi , nrow = 1)
  } else if (!is.matrix(rateParams) ||!is.numeric(rateParams)){
    stop("rateParams must be a numeric matrix")
  } else if (ncol(rateParams)!=nFeaturesRate){
    stop(paste0("ncol(rateParams) must equal the number of coefficients in the formula (",nFeaturesRate,")"))
  } else if(nrow(rateParams) != nEdgeGroup){
    stop("The number of rows in the rateParams matrix must be equal to the number of edgeGroups in the lineageTable")
  }
  if(is.null(piParams)){
    piParams=matrix(1,ncol = nFeaturesPi , nrow = 1)
  } else if (!is.matrix(piParams) ||!is.numeric(piParams) || ncol(piParams)!=nFeaturesPi){
    stop(paste0("piParams must be a numeric matrix with ncol=number of coefficients in the formula (",nFeaturesPi,")"))
  }
  
  ## **End parameter tests** 
  
  ## Get number of sites
  nSites=nrow(covariateTable)
  
  ## Augment piParams with base allele level
  piParamsAug=rbind(matrix(0,ncol=nFeaturesPi),piParams)
  
  ## Create the feature tables
  rateFeatureTable=model.matrix(rateFormula,covariateTable)
  piFeatureTable=model.matrix(rateFormula,covariateTable)
  
  ## Ensure lineageTable edge groups have contiguous integer values starting at 1
  lineageTable[,edgeGroup:=as.integer(as.factor(edgeGroup))]
  
  ## compute pi for all sites
  piAll=exp(piFeatureTable %*%  t(piParamsAug))/rowSums(exp(piFeatureTable %*%  t(piParamsAug)))
  nAlleles=ncol(piAll)  
  ## Create rate parameter matrix where col = branch parameters and row = coefficients
  branchRateParams=apply(tr$edge, 1, function(x) rateParams[lineageTable[child==x[2]]$edgeGroup,])
  ## Compute the rate for all sites
  branchRateAll=exp(rateFeatureTable %*% branchRateParams)

  ## Create matrix to hold simulated data
  simDat=matrix(nrow = nSites,ncol=length(tr$tip.label))
  colnames(simDat)=tr$tip.label
    
  ## Simulate data for each site
  for(i in 1:nSites){
    sitePi=piAll[i,]
    ## Normalize the rate
    temp=matrix(1,ncol=nAlleles,nrow = nAlleles)
    diag(temp)=0
    Q=temp %*% diag(sitePi) ## constuction that guarentees detailed balance: pi_i*q_ij = pi_j*q_ji
    diag(Q)=-rowSums(Q)
    normRate=branchRateAll[i,]/sum(-diag(Q)*sitePi)
    ## Rescale tree edges
    trRescale=tr
    trRescale$edge.length=tr$edge.length*normRate
    ## Simulate data for site i
    simDat[i,] <- ape::rTraitDisc(phy = trRescale,rate=1,k = length(sitePi),freq=sitePi,ancestor = FALSE,
                                  root.value = sample(x=1:length(sitePi),size = 1,prob = sitePi))
  }
  rm(piAll,branchRateParams)
  aData=disCharToProb(simDat,charLevels=1:length(pi))
  return(list(data=aData,nAlleles=nAlleles))
}
