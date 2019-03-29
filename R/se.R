#' Compute standard errors for parameters
#'
#' Estimates the standard error for all non-fixed parameters in the model using the hessian
#' @param model rateModel object
#' @param hess hessian matrix for (not supplied if x is a rateModel object)
#' @name se
#' @return a data.table with the parameter values and standard errors
#' @export
methods::setGeneric("se", function(model, hess = NULL) {
  standardGeneric("se")
})

#' @name se
#' @rdname se
#' @aliases se,rateModel,rateModel-method
methods::setMethod("se", signature(model = "rateModel", hess = "missing"), function(model, hess = NULL) {
  ## Get parameter values
  paramVals <- getParams(model)[which(!model@fixed)]
  ## Compute the hessian
  hess <- numDeriv::hessian(func = phyloGLM:::scaledLL, x = paramVals, model = model, scale = 1)
  out <- se(paramVals, hess)
  return(out)
})

#' @name se
#' @rdname se
#' @aliases se,rateModel,rateModel-method
methods::setMethod("se", signature(model = "rateModel", hess = "matrix"), function(model, hess) {
  if (length(getParams(model)[which(!model@fixed)]) != nrow(hess)) {
    stop("Hessian has more rows than model has unfixed parameters")
  }
  ## Get parameter values
  paramVals <- getParams(model)[which(!model@fixed)]
  ## Invert the hessian
  covM <- solve(-hess)
  stdErr <- sqrt(diag(covM))
  ## Construct the out table to return
  temp <- rbind(
    data.table::data.table(model@phylogeny$getRateIndex(), pType = "rate"),
    data.table::data.table(model@phylogeny$getPiIndex(), pType = "pi")
  )
  ## Adjust idx to idx+1 to make 1-based
  out <- temp[, .(group, name, idx = idx + 1, pType)]
  ## Get params
  out[, value := getParams(model)[idx]]
  data.table::setkey(out, "idx")
  ## Add standard errors
  out[which(!model@fixed), se := stdErr]
  return(out)
})
