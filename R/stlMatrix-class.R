#' Class stlMatrix
#'
#' Class \code{stlMatrix} holds an external pointer to an stl matrix and had various methods
#' to make it easier to access elements
#'
#' @name stlMatrix-class
#' @rdname stlMatrix-class
## #' @include stlMatrixValidityCheck.R
#' @exportClass stlMatrix
methods::setClass("stlMatrix", slots = c(x = "externalptr", nrow = "integer", ncol = "integer"))

#' stlMatrix
#'
#' Contructs an object that holds a pointer to an stl matrix
#' @param x a matrix
#' @name stlMatrix
#' @return an stlMatrix object
stlMatrix <- function(x) {
  if (class(x) != "matrix") {
    stop("x must be a matrix")
  }
  return(new("stlMatrix", x = matrixToStlXptr(x), nrow = nrow(x), ncol = ncol(x)))
}

#' @export
setMethod(
  "[", c("stlMatrix", "numeric", "missing", "ANY"),
  function(x, i, j, ..., drop = TRUE) {
    j <- 0:(x@ncol - 1)
    return(stlMatrixSubset(x@x, i - 1, j))
  }
)

#' @export
setMethod(
  f = "[", signature = c("stlMatrix", "missing", "numeric", "ANY"),
  function(x, i, j, ..., drop = FALSE) {
    i <- 0:(x@nrow - 1)
    return(stlMatrixSubset(x@x, i, j - 1))
  }
)

#' @export
setMethod(
  f = "[", signature = c("stlMatrix", "numeric", "numeric", "ANY"),
  function(x, i, j, ..., drop = FALSE) {
    return(stlMatrixSubset(x@x, i - 1, j - 1))
  }
)
