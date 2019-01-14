#' Discrete characters to probabailites
#'
#' Converts discrete characters labels at each site to probabailites
#' @param x either a matrix or a data.frame of allele values with each column being a species and each row being a site 
#' @param levels a vector of all possible values of the discrete allele label
#' @name disCharToProb
#' @return a list of allele probabilities with one matrix per species
#' @export 
disCharToProb <- function(x,charLevels){
  if(is.null(colnames(x)) || length(unique(colnames(x)))!=ncol(x)){
    colnames(x)=1:ncol(x)
    warning("Columns of x do not have unique column names.")  
  }
  aData=lapply(split(t(x), f =colnames(x)),function(y){
    z=matrix(0,nrow = length(y),ncol=length(charLevels))
    yConv=as.numeric(factor(y,levels = charLevels))
    if(any(is.na(yConv))){
      stop("Invalid character level present in data")
    }
    for(k in 1:length(yConv)){z[k,yConv[k]]=1}
    return(z)
  })
  return(aData)
}
