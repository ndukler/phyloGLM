methods::setGeneric("fit", function(obj,scale=NULL,method=c("l-bfgs-b","mlsl","stogo"),threads=1,control=list()) {
  standardGeneric("fit")
})

#' Fits rate model
#'
#' Fits rate model object and returns fitted model
#' @param obj rateModel
#' @param scale a scale factor to apply to log-likelihood, defaults to -1/nsites
#' @param method Optimization method to use ("l-bfgs-b","mlsl","stogo")
#' @param threads number of threads to use
#' @param control See control from \link[stats]{optim} if using l-bfgs-b, otherwise look at \link[nloptr]{nl.opts}
#' @name fit
#' @return a list incuding information about the optimization, the model object is updated directly (by reference)
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("fit", signature(obj = "rateModel"), function(obj,scale=NULL,method=c("l-bfgs-b","mlsl","stogo"),
                                                                 threads=1,control=list()) {
  ## scale defaults to -1/nsites
  if(is.null(scale)){
    sca=-1/getAlleleData(obj)@data@nrow
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
    cont = list(ndeps=rep(10^-6,sum(!obj@fixed)))
  } else {
    cont = nloptr::nl.opts()
  }
  ## Overwrite defaults when user has supplied values
  for(n in names(control)){
    cont[[n]]=control[[n]]
  } 
  
  if(method=="l-bfgs-b"){
    optMod=optim(par = getParams(obj)[which(!obj@fixed)],fn = scaledLL,obj=obj,scale=sca,method="L-BFGS-B",
               threads=threads,control = cont,hessian = TRUE)
  } else if(method == "mlsl"){
    stop("Unimplemented optimization method specified")
    optMod=nloptr::mlsl(x0=getParams(obj)[which(!obj@fixed)],fn = scaledLL,obj=obj,scale=sca,
                        threads=threads,control = cont)
    optMod$hessian=numDeriv::hessian(func = scaledLL,x=optMod$par,obj=obj,scale=1)
    counts=optMod$iter
  } else if(method == "stogo"){
    stop("Unimplemented optimization method specified")
    optMod=nloptr::stogo(x0=getParams(obj)[which(!obj@fixed)],fn = scaledLL,obj=obj,scale=sca,
                         threads=threads)
    optMod$hessian=numDeriv::hessian(func = scaledLL,x=optMod$par,obj=obj,scale=1)
    counts=optMod$iter
  } else {
    stop("Invalid optimization method specified")
  }
  setParams(obj,optMod$par,which(!obj@fixed)-1)
  return(with(optMod,list(value=value,counts=counts,convergence=convergence,message=message,par=optMod$hessian=hessian)))
})

