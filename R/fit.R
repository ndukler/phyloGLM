methods::setGeneric("fit", function(rateModel,scale=NULL,method=c("l-bfgs-b","mlsl","stogo"),threads=1,control=list(),
                                    log="log.txt") {
  standardGeneric("fit")
})

#' Fits rate model
#'
#' Fits rate model object and returns fitted model
#' @param rateModel rateModel
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
methods::setMethod("fit", signature(rateModel = "rateModel"), function(rateModel,scale=NULL,method=c("l-bfgs-b","mlsl","stogo"),
                                                                 threads=1,control=list(),log="log.txt") {
  ## scale defaults to -1/nsites
  if(is.null(scale)){
    scale=-1/getAlleleData(rateModel)@data@nrow
  }
  ## Check method and set defaults
  if(length(method)>1){
    method="l-bfgs-b"
  } else if(! method %in% c("l-bfgs-b","mlsl","stogo")){
    stop("Invalid optimization method specified")
  }
  
  ## NEED TO FIX SO THAT RATE PARAMETERS CAN BE FIXED INDIVIDUALLY BUT PI IS ALL OR NOTHING PER SITE!!!!!
  ## Where the parameter values are fixed, set lb=ub=value
  ## ub[rateModel@fixed]=x[rateModel@fixed]
  ## lb[rateModel@fixed]=x[rateModel@fixed]
  lb=rep(-10,length(getParams(rateModel)))
  ub=rep(10,length(getParams(rateModel)))
  
  ## Set default control options
  if(method %in% c("l-bfgs-b")){
    cont = list(ndeps=rep(10^-6,sum(!rateModel@fixed)))
  } else {
    cont = nloptr::nl.opts()
  }
  ## Overwrite defaults when user has supplied values
  for(n in names(control)){
    cont[[n]]=control[[n]]
  }
  cont[["trace"]]=1
  
  if(method=="l-bfgs-b"){
    sink(file=log)
    optMod=ucminf::ucminf(par = getParams(rateModel)[which(!rateModel@fixed)],fn = scaledLL,rateModel=rateModel,scale=scale,
                          threads=threads,hessian = TRUE)
    sink()
    optMod$counts=optMod$info[4]
    optMod$hessian=1/scale*optMod$hessian ## revert scaling on hessian
  } else if(method == "mlsl"){
    stop("Unimplemented optimization method specified")
    optMod=nloptr::mlsl(x0=getParams(rateModel)[which(!rateModel@fixed)],fn = phyloGLM:::scaledLL,rateModel=rateModel,scale=scale,
                        threads=threads,control = cont,lower = lb,upper = ub)
    optMod$hessian=numDeriv::hessian(func = scaledLL,x=optMod$par,rateModel=rateModel,scale=1)
    counts=optMod$iter
  } else if(method == "stogo"){
    stop("Unimplemented optimization method specified")
    optMod=nloptr::stogo(x0=getParams(rateModel)[which(!rateModel@fixed)],fn = phyloGLM:::scaledLL,rateModel=rateModel,scale=scale,
                         threads=threads,lower = lb,upper = ub)
    optMod$hessian=numDeriv::hessian(func = scaledLL,x=optMod$par,rateModel=rateModel,scale=1)
    counts=optMod$iter
  } else {
    stop("Invalid optimization method specified")
  }
  setParams(rateModel,optMod$par,which(!rateModel@fixed)-1)
  return(with(optMod,list(value=value,counts=counts,convergence=convergence,message=message,
                          par=optMod$par,hessian=optMod$hessian)))
})

