#' Class rateModel
#'
#' Class \code{rateModel} holds a rate model for a set of alleles
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
#' @name rateModel-class
#' @rdname rateModel-class 
#' @include rateModelValidityCheck.R
#' @include stlMatrix-class.R
#' @importClassesFrom data.table data.table
#' @exportClass rateModel
methods::setClass("rateModel", slots=c(alleleData = "environment",edgeGroups="data.table",
                                       rateFormula="formula",piFormula="formula",rateDM="stlMatrix",
                                       piDM="stlMatrix",phylogeny="ANY",fixed="logical"),
                  validity = rateModelValidityCheck)

#' rateModel
#'
#' Contructs an object that holds a rate model for a set of alleles
#' @param data An alleleData object
#' @param rateFormula A formula that uses the variables in the allele data object to compute turnover rate
#' @param piFormula A formula that uses the variables in the allele data object to compute pi (if NULL uses the same as rateFormula)
#' @param lineageTable A table with columns parent,child,edgeID, and rateGroup, where rateGroup is used to specify how the rates are tied between the branches
#' @name rateModel
#' @return an rateModel object
#' @examples
#'
#' @export
rateModel <- function(data,rateFormula,piFormula=NULL,lineageTable=NULL){
  ## ** Validity checks and default setting** ##
  ## Check that data is of class alleleData
  if(class(data)!="alleleData"){
    stop("data must be of class alleleData")
  }
  ## Check that a rate formula is specified as a formula
  if(class(rateFormula)!="formula"){
    stop("rateFormula must be a formula")
  }
  ## If piFormula is NULL set equal to rateFormula
  if(is.null(piFormula)){
    write("piFormula is not specified, using same formula as rateFormula...")
    piFormula=rateFormula
  }
  ## Check that all coavariates specified in rateFormula are contained in siteInfo
  if(!all(all.vars(rateFormula) %in% colnames(getSiteInfo(data)))){
    stop("Some of the covariates in the rateFormula are not present in the siteInfo of the alleleData object.")
  }
  ## Check that all coavariates specified in piFormula are contained in siteInfo
  if(!all(all.vars(piFormula) %in% colnames(getSiteInfo(data)))){
    stop("Some of the covariates in the piFormula are not present in the siteInfo of the alleleData object.")
  }
  ## Warn user if the formula contains an intercept and provide info on how to avoid 
  if(attr(terms(rateFormula),"intercept")!=0){
    warning("rateFormula has an intercept term. To prevent this specify formula with +0 at the end.")
  }
  ## Warn user if the formula contains an intercept and provide info on how to avoid 
  if(attr(terms(piFormula),"intercept")!=0){
    warning("piFormula has an intercept term. To prevent this specify formula with +0 at the end.")
  }
  ## Checks for lineage table validity
  if(!is.null(lineageTable)){
    if(!is.data.frame(lineageTable)){
      stop("Lineage table must be a data.frame or data.table")
    } else if(!all(c("parent","child","edgeGroup") %in% colnames(lineageTable))){
      stop("Lineage table must contain the columns \'parent\', \'child\', and \'edgeGroup\'")      
    } else if(!setequal(with(getEdgeTable(data),paste0(parent,"-",child)),with(lineageTable,paste0(parent,"-",child)))){
      stop("Edges in the alleleData object and the supplied lineageTable do not match. Run getEdgeTable(data) to view alleleData edges.")  
    } else if(any(table(with(lineageTable,paste0(parent,"-",child)))>1)){
      stop("Duplicated edges in lineageTable")
    }
    ## Ensuring edge groups are integer labeled from 0 to number_of_groups-1
    lineageTable[,edgeGroup:=as.numeric(as.factor(edgeGroup))-1]
  } else {
    ## Create default lineage table
    lineageTable=getEdgeTable(data)
    lineageTable[,edgeGroup:=0]
  }
  
  ## ** Intermediate reformating and computation ** ##
  ## Create environment to hold parameter values and associated indices, etc.
  adEnviron=base::new.env()
  adEnviron$alleleData=data
  
  ## Create the design matrix for the rate variable
  rateDM=stats::model.matrix(rateFormula,getSiteInfo(data))
  ## Create pi design matrix
  if(isTRUE(all.equal.formula(rateFormula,piFormula))){
    piDM=rateDM
  } else {
    piDM=stats::model.matrix(piFormula,getSiteInfo(data))
  }
  
  ## Standardize format for lineageTable and sort by child
  data.table::setcolorder(lineageTable,c("parent","child","edgeGroup"))
  data.table::setkeyv(x = lineageTable,cols = c("child"))
  
  # Create parameter index
  rateP=expand.grid(group=unique(lineageTable$edgeGroup),column=1:ncol(rateDM)-1,
                    stringsAsFactors = FALSE)
  piP=expand.grid(group=2:data@nAlleles-2,column=1:ncol(piDM)-1,stringsAsFactors = FALSE)
  ## Create parameter vector
  params=rep(0,nrow(rateP)+nrow(piP))
  ## Create vector of edge lengths that are indexed by the child id #  
  eL = rep(-1,max(getTree(data)$edge[,2]))
  eL[getTree(data)$edge[,2]] = getTree(data)$edge.length
  ## Get the list of sibblings for every node
  sib=list()
  for(n in unique(as.vector(getTree(data)$edge))){
    sib[[n]]=getSiblings(getTree(data),n)-1
  }
  ## Collect tree info in list
  treeInfo=list(getTree(data)$edge-1, eL,length(getTree(data)$tip.label),sib)
  ## Create edge group vector ordered by child id #
  eG = rep(-1,max(lineageTable$child))
  eG[lineageTable$child]=lineageTable$edgeGroup
  ## Set fixed vector to default to FALSE
  fixed=logical(length(params))
  
  ## ** Object construction ** ##
  methods::new("rateModel",alleleData=adEnviron,edgeGroups=lineageTable,rateFormula=rateFormula,
               piFormula=piFormula,rateDM=stlMatrix(rateDM),piDM=stlMatrix(piDM),
               phylogeny=new(phyloGLM:::phylogeny,params,data.frame(rateP$group,
                    rateP$column,colnames(rateDM)[rateP$column+1]),
                    data.frame(piP$group,piP$col,colnames(piDM)[piP$column+1]),eG,treeInfo),
               fixed=fixed)
}
