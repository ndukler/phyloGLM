#' Plot turnover events on the tree
#'
#' Plots tree with various helpful highlightings
#' @param model rateModel
#' @param aggregate the type of aggregation to perform (either "edge" or "node", see \code{marignalTransitions})
#' @param excludeNodes Indicates nodes that should be excluded from the per-edge aggregated statistic
#' (1-based).
#' @param excludeEdges Indicates edges that should be excluded from the per-node aggregated statistic
#' (edge ids based on integer id of child, use plotTree(model)
#' to view).
#' @param offset How much to offset tip labels by
#' @param xmax Maximum value of x axis in tree plots (Note: if you want different x axis for different tree
#' plots use \link[ggtree]{xlim_expand} on the returned object)
#' @param nodeLabels if TRUE label nodes with numbers
#' @param tipLabels if TRUE label tips with names
#' @name plotElementTransitions
#' @include alleleData-class.R
#' @include rateModel-class.R
#' @rdname plotElementTransitions
#' @examples
#'
#' @export
methods::setGeneric("plotElementTransitions", function(model,aggregate = c("edge","node"),
                                                       excludeNodes = integer(0), excludeEdges = integer(0),
                                                       offset = 0.1, xmax = NULL, nodeLabels = TRUE, tipLabels = TRUE) {
  standardGeneric("plotElementTransitions")
})

#' @name plotElementTransitions
#' @rdname plotElementTransitions
#' @aliases plotElementTransitions,rateModel,rateModel-method
methods::setMethod("plotElementTransitions", signature(model = "rateModel"),
                   function(model, aggregate = c("edge","node"), excludeNodes = integer(0), excludeEdges = integer(0),
                            offset = 0.1, xmax = NULL, nodeLabels = TRUE, tipLabels = TRUE) {
  if(aggregate[1] == "edge"){
    edge_transitions <- marginalTransitions(model, aggregate = "edge",excludeNodes = excludeNodes)
    ## Place data in plottable table
    edge_transition_table <- data.table::rbindlist(lapply(edge_transitions, function(x) {
                                                            data.table::as.data.table(data.table::melt(x))
                                                          }),idcol = "edge")

    data.table::setnames(edge_transition_table, c("Var1", "Var2"), c("ancestral_state", "descendant_state"))
    ## Remove self transitions
    edge_transition_table <- edge_transition_table[ancestral_state != descendant_state,]
    ## Remove edges that should be excluded
    edge_transition_table <- edge_transition_table[!edge %in% excludeEdges,]
    ## Create turnover label
    edge_transition_table[,turnover_label := paste0(ancestral_state,"->",descendant_state)]
    ## Compute percentages for gain/loss events
    total_events = sum(edge_transition_table$value)
    edge_transition_table[,percent_turnovers := value/total_events,by = "turnover_label"]
    ## Reshape edge_transition_table to have one column per turnover type
    turnover_matrix = data.table::dcast(edge_transition_table,edge~turnover_label,value.var = "value")

    ## Merge data to plot gains and losses into tree data
    tree_plot = plotTree(
      model,
      offset = offset,
      nodeLabels = nodeLabels,
      tipLabels = tipLabels,
      xmax = xmax
    )
    tree_plot$data = merge(
      tree_plot$data,
      turnover_matrix,
      by.x = "node",
      by.y = "edge",
      all.x = TRUE
    )
    tree_plot$data$y = tree_plot$data$y * (max(tree_plot$data$x) / max(tree_plot$data$y))

    ## compute properties of the pie charts
    max_rad = 0.15
    subset_columns = grep("node|x|y|->|branch.length", colnames(tree_plot$data))
    pi_info = as.data.table(tree_plot$data[, subset_columns])
    pi_info[, x := (2 * x - branch.length) / 2] ## update x so it is the midpoint of the branch
    transition_columns = grep("->", colnames(pi_info))
    total_turnovers = sum(pi_info[, ..transition_columns], na.rm = TRUE)
    edge_turnovers = rowSums(pi_info[, ..transition_columns])
    pi_info[, area := ..edge_turnovers / ..total_turnovers]
    pi_info[, radius := sqrt(area / 3.1415)] ## compute radius from area
    pi_info[, scaled.radius := radius / max(radius, na.rm = TRUE) * max_rad] ## scale radius based on max radius

    pi_info = pi_info[!is.na(radius), ]
    turnovers = scales::comma(round(colSums(pi_info[, ..transition_columns], na.rm =
                                              TRUE))) ## total promoter gains


    final_plot = tree_plot +
      scatterpie::geom_scatterpie(
        mapping = ggplot2::aes(x = x, y = y, r = scaled.radius),
        data = pi_info,
        cols = colnames(pi_info)[transition_columns]
      ) +
      ggplot2::scale_fill_discrete(name = "Turnover type") +
      ggplot2::scale_x_continuous(
        name = "Time",
        labels = function(x)
          round(abs(x - max(tree_plot$data$x)), 2)
      ) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 16),
        legend.position = "right",
        legend.background = ggplot2::element_blank()
      ) +
      ggplot2::guides(colour = FALSE) +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 15),
        axis.title = ggplot2::element_text(size = 16)
      ) +
      scatterpie::geom_scatterpie_legend(
        radius = pi_info$scaled.radius[pi_info$scaled.radius > 0],
        x = 0,
        y = max(pi_info$y) - max(pi_info$scaled.radius),
        labeller = function(x)
          scales::percent((x / max_rad * max(pi_info$radius)) ^ 2 * 3.1415),
        n = 4
      )


  } else if (aggregate[1] == "node"){
    stop("Per-node element transition plotting not implemented")
  } else {
    stop("Invalid type of aggregation specified. Valid types are \'edge\' and \'node\'.")
  }
  return(final_plot)
})