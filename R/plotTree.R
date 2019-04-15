#' Tree visualization
#'
#' Plots tree with various helpful highlightings
#' @param model alleleData
#' @param offset How much to offset tip labels by
#' @param xmax Maximum value of x axis in tree plots (Note: if you want different x axis for different tree 
#' plots use \link[ggtree]{xlim_expand} on the returned object)
#' @param nodeLabels if TRUE label nodes with numbers
#' @param nodeLabels if TRUE label tips with names
#' @name plotTree
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @rdname plotTree
#' @examples
#' 
#' @export
methods::setGeneric("plotTree", function(model,...) {
  standardGeneric("plotTree")
})

#' @name plotTree
#' @rdname plotTree
#' @aliases plotTree,rateModel,rateModel-method
methods::setMethod("plotTree", signature(model = "rateModel"), function(model,offset=0.1,xmax=NULL,nodeLabels=TRUE,tipLabels=TRUE) {
  tempTree=getTree(model)
  g <- ggtree::ggtree(tempTree)+
    ggtree::theme_tree2()
  if(tipLabels){g <- g+ggtree::geom_tiplab(show.legend=FALSE,color="black",offset=offset)}
  if(nodeLabels){g <- g+ggplot2::geom_label(ggplot2::aes(label=node), hjust=0.5)}
  ## Copied over code from ggtree since xlim_expand not available for older versions of R/Bioconductor
  if(!is.null(xmax)) {
    dummy <- data.frame(x=xmax, .panel='Tree')
    g=g+geom_blank(aes_(x=~x), dummy, inherit.aes = FALSE)
  }
  ## Add color based on edge group
  temp <- data.table::copy(model@edgeGroups)
  setkey(temp,"child")
  g$data$group = temp[.(g$data$node),]$edgeGroup
  if (requireNamespace("randomcoloR", quietly = TRUE)){
    g=g+ggplot2::aes(color=index)+
      ggplot2::scale_color_manual(values=randomcoloR::distinctColorPalette(length(levels(g$data$group))),breaks=levels(temp$group))
  } else {
    g=g+ggplot2::aes(color=group)
  }
  return(g)
})
