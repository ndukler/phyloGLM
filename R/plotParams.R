#' Plot parameter values
#'
#' Plots parameter values and standard errors
#' @param obj rateModel
#' @name plotParams
#' @return a data.table with the parameter values and standard errors
#' @include rateModel-class.R
#' @export 
methods::setGeneric("plotParams", function(obj,type=c("bar","dot")) {
  standardGeneric("plotParams")
})

#' @name plotParams
#' @rdname plotParams
#' @aliases plotParams,rateModel,rateModel-method 
methods::setMethod("plotParams", signature(obj = "rateModel"), function(obj,type="dot") {
  # ## Check one type is selected
  # if(type == c("bar","dot")){
  #   type="bar"
  # } else if(length(type)!=1){
  #   stop("Please specify a single plot type option")
  # } else if (!type %in% c("bar","dot") ){
  #   stop("Invalid plot type selected")
  # }
  # 
  ## Get parameter values and standard errors
  seTab=se(obj)
  ## Create bar plot for parameters
  if("dot"){
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
                                                ymin = value - se,
                                                ymax = value + se, 
                                                color= factor(group+2))) +
      ggplot2::geom_point(size=2, position=pos, stat="identity") +
      ggplot2::geom_errorbar(position=pos, ggplot2::aes(width=0.2))+
      cowplot::theme_cowplot()+
      ggplot2::guides(color=ggplot2::guide_legend(title="Allele"))+
      ggplot2::xlab("Parameter")+
      ggplot2::ylab("Value")
    gAll <- cowplot::plot_grid(rateG,piG,ncol=1)
  }
  return(gAll)
})