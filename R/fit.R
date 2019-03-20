methods::setGeneric("fit", function(model,scale=NULL,method=c("l-bfgs-b","mlsl","stogo"),threads=1,control=list(),
                                    log="log.txt") {
  standardGeneric("fit")
})

#' Fits rate model
#'
#' Fits rate model object and returns fitted model
#' @param model rateModel
#' @param scale a scale factor to apply to log-likelihood, defaults to -1/nsites
#' @param method Optimization method to use ("bfgs","mlsl","stogo")
#' @param threads number of threads to use
#' @param control See control from \link[stats]{optim} if using l-bfgs-b, otherwise look at \link[nloptr]{nl.opts}
#' @name fit
#' @return a list incuding information about the optimization, the model object is updated directly (by reference)
#' @include rateModel-class.R
#' @examples
#' 
#' @export
methods::setMethod("fit", signature(model = "rateModel"), function(model,scale=-1,method=c("bfgs","mlsl","stogo"),
                                                                 threads=1,control=list(),log="log.txt") {
  ## scale defaults to -1/nsites
  if(!is.numeric(scale)){
    stop("scale must be a numeric value")
  }
  ## Check method and set defaults
  if(length(method)>1){
    method="bfgs"
  } else if(! method %in% c("bfgs","mlsl","stogo")){
    stop("Invalid optimization method specified")
  }
  
  ## NEED TO FIX SO THAT RATE PARAMETERS CAN BE FIXED INDIVIDUALLY BUT PI IS ALL OR NOTHING PER SITE!!!!!
  ## Where the parameter values are fixed, set lb=ub=value
  ## ub[model@fixed]=x[model@fixed]
  ## lb[model@fixed]=x[model@fixed]
  lb=rep(-5,length(getParams(model)))
  ub=rep(5,length(getParams(model)))
  
  ## Set default control options
  if(method %in% c("bfgs")){
    cont = list(ndeps=rep(10^-6,sum(!model@fixed)))
  } else {
    cont = nloptr::nl.opts()
  }
  ## Overwrite defaults when user has supplied values
  for(n in names(control)){
    cont[[n]]=control[[n]]
  }
  cont[["trace"]]=1
  
  if(method=="bfgs"){
    sink(file=log)
    optMod=optim(par = getParams(model)[which(!model@fixed)],fn = scaledLL,gr = phyloGrad,model=model,scale=scale,
                          threads=threads,method="BFGS",hessian = FALSE)
    sink()
    setParams(model,optMod$par,which(!model@fixed)-1)
    # optMod$counts=optMod$iter
    optMod$hessian=numDeriv::hessian(func = ll,x=optMod$par,model=model)
  } else if(method == "mlsl"){
    stop("Unimplemented optimization method specified")
    optMod=nloptr::mlsl(x0=getParams(model)[which(!model@fixed)],fn = phyloGLM:::scaledLL,model=model,scale=scale,
                        threads=threads,control = cont,lower = lb,upper = ub)
    setParams(model,optMod$par,which(!model@fixed)-1)
    optMod$hessian=numDeriv::hessian(func = scaledLL,x=optMod$par,model=model,scale=1)
    optMod$counts=optMod$iter
  } else if(method == "stogo"){
    optMod=nloptr::stogo(x0=getParams(model)[which(!model@fixed)],fn = phyloGLM:::scaledLL,model=model,scale=scale,
                         threads=threads,lower = lb,upper = ub)
    setParams(model,optMod$par,which(!model@fixed)-1)
    optMod$hessian=numDeriv::hessian(func = scaledLL,x=optMod$par,model=model,scale=1)
    counts=optMod$iter
  } else {
    stop("Invalid optimization method specified")
  }
  return(with(optMod,list(value=value,counts=counts,convergence=convergence,message=message,
                          par=optMod$par,hessian=optMod$hessian)))
})

