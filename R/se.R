#' Compute standard errors for parameters
#'
#' Estimates the standard error for all non-fixed parameters in the model using the hessian
#' @param x rateModel object
#' @param hess hessian matrix for (not supplied if x is a rateModel object)
#' @name se
#' @return a data.table with the parameter values and standard errors
#' @export 
methods::setGeneric("se", function(obj,hess=NULL) {
  standardGeneric("se")
})

#' @name se
#' @rdname se
#' @aliases se,rateModel,rateModel-method 
methods::setMethod("se", signature(obj = "rateModel",hess="missing"), function(obj,hess=NULL) {
  ## Get parameter values
  paramVals=getParams(obj)[which(!obj@fixed)]
  ## Compute the hessian
  hess=numDeriv::hessian(func = phyloGLM:::scaledLL,x=paramVals,obj=obj,scale=1)
  out=se(paramVals,hess)
  return(out)
})

#' @name se
#' @rdname se
#' @aliases se,rateModel,rateModel-method 
methods::setMethod("se", signature(obj = "rateModel",hess="matrix"), function(obj,hess) {
  ## Get parameter values
  paramVals=getParams(obj)[which(!obj@fixed)]
  ## Invert the hessian
  covM<-solve(-hess)
  stdErr<-sqrt(diag(covM))
  ## Construct the out table to return
  temp=rbind(data.table::data.table(obj@phylogeny$getRateIndex(),pType="rate"),
             data.table::data.table(obj@phylogeny$getPiIndex(),pType="pi"))
  ## Adjust idx to idx+1 to make 1-based
  out=temp[,.(group,name,idx=idx+1,pType)]
  ## Get params
  out[,value:=getParams(obj)[idx]]
  data.table::setkey(out,"idx")
  ## Add standard errors
  out[which(!obj@fixed),se:=stdErr]
  return(out)
})