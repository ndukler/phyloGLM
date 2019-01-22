#' Compute standard errors for parameters
#'
#' Estimates the standard error for all non-fixed parameters in the model using the hessian
#' @param obj rateModel
#' @name se
#' @return a data.table with the parameter values and standard errors
#' @include rateModel-class.R
#' @export 
methods::setGeneric("se", function(obj) {
  standardGeneric("se")
})

#' @name se
#' @rdname se
#' @aliases se,rateModel,rateModel-method 
methods::setMethod("se", signature(obj = "rateModel"), function(obj) {
  ## Get parameter values
  paramVals=getParams(obj)[which(!obj@fixed)]
  ## Compute the hessian
  hess=numDeriv::hessian(func = phyloGLM:::scaledLL,x=paramVals,obj=obj,scale=1)
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