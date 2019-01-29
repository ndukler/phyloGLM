#' Simulate data
#'
#' Simulates data without an error model
#' @param tr tree
#' @param covariateTable
#' @param rateFormula
#' @param rateParams
#' @param piFormula
#' @param piParams
#' @param lineageTable
#' @name simulateSites
#' @return a list
#' @export

simulateSites <- function(tr,covariateTable,rateFormula,rateParams,piFormula=NULL,piParams=NULL,lineageTable=NULL){
  ## A bunch of checks 
  
  
  ## Normalize the stationary frequency
  pi=pi/sum(pi)
  nAlleles=length(pi)
  ## Normalize the rate
  temp=matrix(1,ncol=nAlleles,nrow = nAlleles)
  diag(temp)=0
  Q=temp %*% diag(pi) ## constuction that guarentees detailed balance: pi_i*q_ij = pi_j*q_ji
  diag(Q)=-rowSums(Q)
  normRate=rate/sum(-diag(Q)*pi)
  ## Rescale tree edges
  trRescale=tr
  trRescale$edge.length=tr$edge.length*normRate
  ## Create matrix to hold simulated data
  simDat=matrix(nrow = nSites,ncol=length(tr$tip.label))
  colnames(simDat)=tr$tip.label
  for(i in 1:nSites){
    simDat[i,] <- ape::rTraitDisc(phy = trRescale,rate=1,k = length(pi),freq=pi,ancestor = FALSE,
                                  root.value = sample(x=1:length(pi),size = 1,prob = pi))
  }
  aData=disCharToProb(simDat,charLevels=1:length(pi))
  ## Create edgeGroup table
  eTab=data.table::data.table(parent=trRescale$edge[,1],child=trRescale$edge[,2],edgeID=paste(trRescale$edge[,1],trRescale$edge[,2],sep="-"),
                              edgeGroup=paste0("e",as.numeric(factor(rate))))
  ## Map of edge groups to rates
  rateMap=unique(data.table::data.table(edgeGroup=paste0("e",as.numeric(factor(rate))),rate))
  return(list(data=aData,tr=tr,edgeTable=eTab,rateMap=rateMap,pi=pi))
}