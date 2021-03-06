% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal.R
\name{marginalTransitions}
\alias{marginalTransitions}
\title{Expected transitions per site}
\arguments{
\item{model}{rateModel object}

\item{aggregate}{How expected numbers of transitions should be aggregated (within edge or node)}

\item{excludeNodes}{Indicates nodes that should be excluded from the per-edge aggregated statistic 
(1-based). Only usable with aggregate = "edge".}

\item{excludeEdges}{Indicates edges that should be excluded from the per-node aggregated statistic 
Only usable with aggregate = "node". (edge ids based on integer id of child, use plotTree(model) 
to view).}

\item{threads}{number of threads to use}
}
\value{
a list of per node/edge matricies containing the expected number of allelic transitions. The
rows are the ancestral allele and the columns are the descendant allele. The list ids are a character
string that are either a site or edge (child) id depending on how you aggregated
}
\description{
Compute the expected number of transitions aggregated across either nodes or edges
}
