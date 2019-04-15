#' Plot parameter values
#'
#' Plots parameter values and standard errors
#' @param x rateModel object or the data.frame from se(rateModel)
#' @param hess the hessian for the model (optional)
#' @param alpha plot the confidence interval of the parameters from [alpha/2,1-alpha/2]
#' @name plotParams
#' @return a list of two plots, one for the rate parameters, one for the pi parameters
#' @include rateModel-class.R
#' @export
methods::setGeneric("plotParams", function(x, hess = NULL, alpha = 0.05, subset = NULL) {
  standardGeneric("plotParams")
})

#' @name plotParams
#' @rdname plotParams
#' @aliases plotParams,rateModel,rateModel-method
methods::setMethod("plotParams", signature(x = "rateModel"), function(x, hess = NULL, alpha = 0.05, subset=NULL) {
  if (!is.null(hess) && !is.matrix(hess)) {
    stop("Hessian must be a matrix")
  }
  ## Get parameter values and standard errors
  seTab <- se(model = x, hess = hess)
  ## Create bar plot for parameters
  gAll <- plotParams(seTab, hess = NULL, alpha = alpha, subset = subset)
  return(gAll)
})

#' @name plotParams
#' @rdname plotParams
#' @aliases plotParams,data.frame,data.frame-method
methods::setMethod("plotParams", signature(x = "data.frame"), function(x, hess = NULL, alpha = 0.05, subset = NULL) {
  ## Check if the data.frame has the correct columns
  if (!all(c("group", "name", "idx", "pType", "value", "se") %in% colnames(x))) {
    stop("Incorrectly formated parameter table. Try running se() on the rateModel object to get the correct table.")
  }
  ## Compute factor to multiply se by to compute ci
  z_star <- qnorm(p = c(alpha / 2, 1 - alpha / 2))
  ## If subset isn't NULL, subset the variables that will be plotted
  if(!is.null(subset)){
    if(any(!subset %in% x$name)){
      not_incl = subset[!subset %in% x$name]
      warning(paste("Subset elements not present in model:",paste(not_incl, collapse = ",")))
    }
    x=x[name %in% subset,]
    if(nrow(x)==0){
      stop("No parameters remaining after subset")
    }
  }
  
  ## Create bar plot for parameters
  pos <- ggplot2::position_dodge(width = 0.25)
  rateG <- ggplot2::ggplot(x[pType == "rate"], ggplot2::aes(
    x = name,
    y = value,
    ymin = value + se * z_star[1],
    ymax = value + se * z_star[2],
    color = factor(group)
  )) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_point(size = 2, position = pos, stat = "identity", shape = 18, size = 1) +
    ggplot2::geom_errorbar(position = pos, ggplot2::aes(width = 0.2)) +
    cowplot::theme_cowplot() +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Edge Group")) +
    ggplot2::xlab("Parameter") +
    ggplot2::ylab("Value") +
    ggplot2::ggtitle("Rate parameters") +
    ggplot2::theme(panel.grid.major = ggplot2::element_line("black", 0.1, "dashed"))
  piG <- ggplot2::ggplot(x[pType == "pi"], ggplot2::aes(
    x = name,
    y = value,
    ymin = value + se * z_star[1],
    ymax = value + se * z_star[2],
    color = factor(group + 2)
  )) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_point(size = 2, position = pos, stat = "identity", shape = 18, size = 1) +
    ggplot2::geom_errorbar(position = pos, ggplot2::aes(width = 0.2)) +
    cowplot::theme_cowplot() +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Allele")) +
    ggplot2::xlab("Parameter") +
    ggplot2::ylab("Value") +
    ggplot2::ggtitle("Allele stationary distribution parameters") +
    ggplot2::theme(panel.grid.major = ggplot2::element_line("black", 0.1, "dashed"))
  return(list(rate = rateG, pi = piG))
})
