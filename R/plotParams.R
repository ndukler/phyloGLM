#' Plot parameter values
#'
#' Plots parameter values and standard errors
#' @param x rateModel object or the data.frame from se(rateModel)
#' @name plotParams
#' @return a data.table with the parameter values and standard errors
#' @include rateModel-class.R
#' @export 
methods::setGeneric("plotParams", function(x) {
  standardGeneric("plotParams")
})

#' @name plotParams
#' @rdname plotParams
#' @aliases plotParams,rateModel,rateModel-method 
methods::setMethod("plotParams", signature(x = "rateModel"), function(x) {
  ## Get parameter values and standard errors
  seTab=se(obj)
  ## Create bar plot for parameters
  gAll <- plotParams(seTab) 
  return(gAll)
})

#' @name plotParams
#' @rdname plotParams
#' @aliases plotParams,data.frame,data.frame-method 
methods::setMethod("plotParams", signature(x = "data.frame"), function(x) {
  ## Check if the data.frame has the correct columns
  if(!all(c("group","name","idx","pType","value","se") %in% colnames(x))){
    stop("Incorrectly formated parameter table. Try running se() on the rateModel object to get the correct table.")
  }
  ## Create bar plot for parameters
  pos = ggplot2::position_dodge(width=0.25)
  rateG = ggplot2::ggplot(seTab[pType=="rate"], ggplot2::aes(x=name,
                                                             y=value,
                                                             ymin = value - se,
                                                             ymax = value + se, 
                                                             color= factor(group))) +
    ggplot2::geom_point(size=2, position=pos, stat="identity") +
    ggplot2::geom_errorbar(position=pos, ggplot2::aes(width=0.2))+
    cowplot::theme_cowplot()+
    ggplot2::guides(color=ggplot2::guide_legend(title="Edge Group"))+
    ggplot2::xlab("Parameter")+
    ggplot2::ylab("Value")
  piG = ggplot2::ggplot(seTab[pType=="pi"], ggplot2::aes(x=name,
                                                         y=value,
                                                         ymin = value - 1.96*se,
                                                         ymax = value + 1.96*se, 
                                                         color= factor(group+2))) +
    ggplot2::geom_point(size=2, position=pos, stat="identity") +
    ggplot2::geom_errorbar(position=pos, ggplot2::aes(width=0.2))+
    cowplot::theme_cowplot()+
    ggplot2::guides(color=ggplot2::guide_legend(title="Allele"))+
    ggplot2::xlab("Parameter")+
    ggplot2::ylab("Value")
  gAll <- cowplot::plot_grid(rateG,piG,ncol=1)
  return(gAll)
})