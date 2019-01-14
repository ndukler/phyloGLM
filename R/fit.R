methods::setGeneric("fit", function(obj,scale=NULL,method=c("l-bfgs-b","mlsl","stogo"),control=list()) {
  standardGeneric("fit")
})

#' Fits rate model
#'
#' Fits rate model object and returns fitted model
#' @param obj rateModel
#' @param scale a scale factor to apply to log-likelihood, defaults to -1/nsites
#' @param method Optimization method to use ("l-bfgs-b","mlsl","stogo")
#' @param control See control from \link[stats]{optim} if using l-bfgs-b, otherwise look at \link[nloptr]{nl.opts}
#' @name fit
#' @return a list incuding a rateModel object and information about the optimization
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("fit", signature(obj = "rateModel"), function(obj,scale=NULL,method=c("l-bfgs-b","mlsl","stogo"),control=list()) {
  ## scale defaults to -1/nsites
  if(is.null(scale)){
    sca=-1/nrow(getAlleleData(obj)@data)
  }
  ## Check method and set defaults
  if(length(method)>1){
    method="l-bfgs-b"
  } else if(! method %in% c("l-bfgs-b","mlsl","stogo")){
    stop("Invalid optimization method specified")
  }
  
  ## NEED TO FIX SO THAT RATE PARAMETERS CAN BE FIXED INDIVIDUALLY BUT PI IS ALL OR NOTHING PER SITE!!!!!
  ## Where the parameter values are fixed, set lb=ub=value
  ## ub[obj@fixed]=x[obj@fixed]
  ## lb[obj@fixed]=x[obj@fixed]
  
  ## Set default control options
  if(method %in% c("l-bfgs-b")){
    cont = list(ndeps=rep(10^-6,length(y)))
  } else {
    cont = nloptr::nl.opts()
  }
  ## Overwrite defaults when user has supplied values
  for(n in names(control)){
    cont[[n]]=control[[n]]
  } 
  
  if(method=="l-bfgs-b"){
    optMod=optim(y,fn = epiAllele:::scaledLL,obj=obj,scale=sca,stickParams=TRUE,lower = lb,upper = ub,method="L-BFGS-B",
               control = cont)
  } else if(method == "mlsl"){
    optMod=nloptr::mlsl(x0=y,fn = epiAllele:::scaledLL,lower = lb,upper = ub,obj=obj,scale=sca,stickParams=TRUE,
                        control = cont)
    counts=optMod$iter
  } else if(method == "stogo"){
    optMod=nloptr::stogo(x0=y,fn = epiAllele:::scaledLL,lower = lb,upper = ub,obj=obj,scale=sca,stickParams=TRUE,
                         control = cont)
    counts=optMod$iter
  } else {
    stop("Invalid optimization method specified")
  }
  probs=multiStickToProb(optMod$par[(rateEnd+1):length(optMod$par)],width=getAlleleData(obj)@nAlleles-1)
  setParamValue(obj2,i = 1:length(obj@params),value = c(optMod$par[1:rateEnd],probs))
  return(with(optMod,list(model=obj2,value=value,counts=counts,convergence=convergence,message=message)))
})

