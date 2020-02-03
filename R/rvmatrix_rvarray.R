
#' @rdname rvarray
#' @param nrow the desired number of rows.
#' @param ncol the desired number of columns.
#' @param byrow logical. If \code{FALSE} (the default) the matrix is filled by
#' columns, otherwise the matrix is filled by rows.
#' @export
rvmatrix <- function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL) {
  data <- as.vector(data)
  if (missing(nrow)) 
      nrow <- ceiling(length(data)/ncol)
  else if (missing(ncol)) 
      ncol <- ceiling(length(data)/nrow)
  if (byrow) {
    X <- t(rvarray(data, c(ncol, nrow), dimnames=dimnames))
  } else {
    X <- rvarray(data, c(nrow, ncol), dimnames=dimnames)
  }
  return(X)
}



#' Matrices and Arrays of Random Vectors
#' 
#' Arrange a given random vector into a matrix or array form.
#' 
#' These are 'rv' compatible versions of the functions \code{\link{matrix}} and
#' \code{\link{array}}.
#' 
#' The function \code{rvmatrix} generates the random variable matrix via an
#' \code{rvarray} call.
#' 
#' The \code{rvarray} function calls first \code{array} to set the dimensions
#' of the argument \code{data} and then coerces the resulting array object to
#' an 'rv' object.
#' 
#' @aliases rvmatrix rvarray is.matrix.rv as.matrix.rv
#' @param data an optional data vector.
#' @param dimnames A dimnames attribute for the matrix: a list of length 2
#' giving the row and column names respectively.
#' @param dim the dim attribute for the array to be created, that is a vector
#' of length one or more giving the maximal indices in each dimension.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso To plot random matrices, see \code{\link{mlplot}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   n.rows <- 3; n.cols <- 4; n <- (n.rows*n.cols)
#'   mu.true <- rnorm(1:n.rows, mean=1:n.rows, sd=1)
#'   theta <- rvmatrix(rvnorm(n=n.cols, mean=mu.true, sd=0.5), nrow=n.rows)
#'   col.labels <- paste("Time", 1:n.cols, sep=":")
#'   row.labels <- paste("Unit", 1:n.rows, sep=":")
#'   dimnames(theta) <- list(row.labels, col.labels)
#'   print(theta)
#'   print(E(theta))  
#' 
#' @export rvarray
rvarray <- function (data = NA, dim = length(data), dimnames = NULL) {
  as.rv(array(data = data, dim = dim, dimnames = dimnames))
}

#' @noRd
#' @method is.matrix rv
#' @param x an R object.
is.matrix.rv <- function (x) {
  dx <- dim(x)
  return((!is.null(dx)) && length(dx)==2)
}

as.matrix.rv <- function (x, ...) {
  if (is.matrix(x)) {
    return(x)
  }
  dn <- if (!is.null(names(x))) list(names(x), NULL) else NULL
  rvarray(x, dim=c(length(x), 1), dimnames=dn)
}

