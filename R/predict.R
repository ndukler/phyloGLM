#' Predict phylogenetic parameter values
#'
#' Predicts the rate and stationary distribution per site
#' @param object rateModel object
#' @param data a new set of genomic covariate values (siteInfo data)
#' @param rateDM allows direct setting of the design matrix for rate covariates 
#' (will override data, reccomended for advanced users only)
#' @param piDM allows direct setting of the design matrix for pi covariates 
#' (will override data, reccomended for advanced users only)
#' @return a list containing both the rates and the stationary distribution at each site
#' @name predict
#' @include rateModel-class.R
#' @rdname predict
#' @export

methods::setMethod("predict", signature(object = "rateModel"), function(object,data = NULL,rateDM = NULL, 
                                                                        piDM = NULL,...) {
  ######################################
  ## Option handling and error checking
  ######################################
  rateDM_new <- NULL
  piDM_new <- NULL
  current.na.action <- options('na.action')
  options(na.action = 'na.pass')
  ## Check that either is.null(data) or is.null(DM)
  if(!is.null(data)){
    if(!is.data.frame(data)){
      stop("data must be a data.frame")
    } else {
      ## Ensure all strings coverted to factors
      data <- data.table::data.table(data, stringsAsFactors = T)
      ## Match types (and where appropriate levels) of matched columns to old model
      for(i in intersect(colnames(data),colnames(getSiteInfo(object)))){
        ## Check that types match
        if(class(data[[i]])!=class(getSiteInfo(object)[[i]])){
          if(is.numeric(data[[i]]) & is.numeric(getSiteInfo(object)[[i]])){
            warning(paste0("Compatible but different data types of column (",i,"): prior data is (",
                           class(getSiteInfo(object)[[i]]),") and new data (",class(data[[i]]),")"))
          } else {
            stop(paste0("Incompatible data types of column ",i,": prior data is (",
                     class(getSiteInfo(object)[[i]]),") and new data (",class(data[[i]]),")"))
          }
        }
        ## If factor, match levels
        if(class(data[[i]])=="factor"){
          data[[i]] = factor(as.character(data[[i]]),levels=levels(getSiteInfo(object)[[i]]))
          if(any(is.na(data[[i]]))){
            stop(paste("NA present in",i,", after factor level coercion to match levels of siteInfo for object."))  
          }
        }
      }
      ## Expand the new data out into design matricies
      rateDM_new <- stats::model.matrix(object@rateFormula, data)
      piDM_new <- stats::model.matrix(object@piFormula, data)
    }
  } 
  ## Overwrite design matricies if desired
  if (!is.null(rateDM)){
    rateDM_new <- rateDM
  }  
  if (!is.null(piDM)){
    piDM_new <- piDM
  }
  options('na.action' = current.na.action$na.action)
  
  ## If rateDM or piDM remain NULL throw error
  if(is.null(rateDM_new)){
    stop("Insufficient information given to construct rateDM. Either provide data or rateDM.")
  }
  if(is.null(piDM_new)){
    stop("Insufficient information given to construct piDM. Either provide data or piDM.")
  }
  
  ## Check that all elements in the design matrix are finite
  if (any(!is.finite(rateDM_new))) {
    stop("Non-numeric/Non-finite elements in the rate design matrix")
  }
  if (any(!is.finite(piDM_new))) {
    stop("Non-numeric/Non-finite elements in the pi design matrix")
  }
  
  ## Check that rate and piDM have the same number of sites
  if(nrow(rateDM_new)!=nrow(piDM_new)){
    stop("Different numbers of sites in rate and design matricies")  
  }
  
  ## Check that rate and piDM have the same number of columns as the design matricies 
  ## in the current model
  if(ncol(rateDM_new)!=object@rateDM@ncol){
    stop("Mismatch in number of columns for previous and new design matricies. 
         May be due to a mismatch in factor levels.")  
  }
  
  ## Convert to stl matricies
  rateDM_ptr <- phyloGLM:::stlMatrix(rateDM_new)
  piDM_ptr <- phyloGLM:::stlMatrix(piDM_new)
  
  ######################################
  ## Predict rates and pi
  ######################################
  sites <- 1:rateDM_ptr@nrow
  ## Select a representitive child for each edge group and adjust index
  repChildren <- object@edgeGroups[!duplicated(object@edgeGroups$edgeGroup), .(child, edgeGroup)]
  ## compute rates
  out_rate <- data.table::as.data.table(expand.grid(edgeGroup = as.integer(repChildren$edgeGroup), site = as.integer(sites)))
  out_rate <- merge(repChildren, out_rate, by = "edgeGroup")
  out_rate[, rate := object@phylogeny$rateV(child - 1, sites - 1, rateDM_ptr@x), by = c("child")]
  out_rate[, child := NULL]
  ## Compute stationary distribution
  out_pi <- data.table::rbindlist(lapply(as.list(sites), function(x) data.table::data.table(
    site = x, allele = 1:getAlleleData(object)@nAlleles,
    object@phylogeny$pi(piDM_ptr[x, ])
  )))
  ## Return predicted values in list
  return(list(rate=out_rate,pi=out_pi))
})
