#' Read in *.paf file
#'
#' @param f path to *.paf file
#' @param subset vector of site labels to be included (defaults to all sites)
#' @name readPAF
#' @return a list of allele probability matricies and a list of the site labels for each row of the matrix
#' @export
readPAF <- function(f,subset){
  ## Check that file exists
  if(!file.exists(f)){
    stop(paste("File",f,"does not exist"))
  }
  ## Check that format has at least 4 columns with fist two named "species" and "site"
  d=data.table::fread(f,nrows = 0)
  if(ncol(d) < 4){
    stop("PAF files must have at least 4 columns, \'site\', \'species\', then at least two alleles.")
  }
  if(any(colnames(d)[1:2]!=c("site","species"))){
    stop("The first two columns of a PAF file must be \'site\' and \'species\' respectively.")
  }
  ## Read in file
  d=data.table::fread(f,key=c("site"))
  ## Check that all sites are present for all species
  allSites=unique(d$site)
  dCheck=d[,setequal(site,allSites),by="species"]
  if(any(!dCheck$V1)){
    stop(paste("Species", paste(dCheck[V1==FALSE]$species,collapse = ","),"are missing sites."))
  }
  ## Subset sites if needed
  if(!is.null(subset)){
    d=d[subset]
  }
  finalSites=unique(d$site)
  ## Convert to species indexed list
  dList=split(d,by="species",keep.by = FALSE)
  dList=lapply(dList,function(x) {data.table::setkey(x,"site");as.matrix(x[finalSites,-1])})
  return(list(dataList=dList,sites=finalSites))
}