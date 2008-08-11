\name{rvarray}
\title{Matrices and Arrays of Random Vectors}
\alias{rvmatrix}
\alias{rvarray}
\alias{is.matrix.rv}
\alias{as.matrix.rv}
%\alias{is.array.rv}
%\alias{as.array.rv}
\usage{
 rvarray(data = NA, dim = length(data), dimnames = NULL) 
 rvmatrix(data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL)
 \method{is.matrix}{rv}(x)
 \method{as.matrix}{rv}(x, \dots)
% \method{is.array}{rv}(x)
% \method{as.array}{rv}(x)
}
\description{
  Arrange a given random vector into a matrix or array form.

  These are 'rv' compatible versions of the functions 
  \code{\link{matrix}} and  \code{\link{array}}.
}
\arguments{
  \item{data}{an optional data vector.}
  \item{nrow}{the desired number of rows.}
  \item{ncol}{the desired number of columns.}
  \item{byrow}{logical. If \code{FALSE} (the default) the matrix is filled by columns, otherwise the matrix is filled by rows.}
  \item{dimnames}{A dimnames attribute for the matrix: a list of length 2 giving the row and column names respectively.}
  \item{dim}{the dim attribute for the array to be created, that is a vector of length one or more giving the maximal indices in each dimension.}
  \item{\dots}{arguments passed to other methods}
  \item{x}{an R object.}
}
\details{
  The function \code{rvmatrix} generates the random variable matrix via an
  \code{rvarray} call.

  The \code{rvarray} function calls first \code{array} to set the dimensions
  of the argument \code{data} and then coerces the resulting array object to an 'rv' object.
}
\seealso{
  To plot random matrices, see \code{\link{mlplot}}.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  n.rows <- 3; n.cols <- 4; n <- (n.rows*n.cols)
  mu.true <- rnorm(1:n.rows, mean=1:n.rows, sd=1)
  theta <- rvmatrix(rvnorm(n=n.cols, mean=mu.true, sd=0.5), nrow=n.rows)
  col.labels <- paste("Time", 1:n.cols, sep=":")
  row.labels <- paste("Unit", 1:n.rows, sep=":")
  dimnames(theta) <- list(row.labels, col.labels)
  print(theta)
  print(E(theta))  
}

\keyword{classes}

