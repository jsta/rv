\name{matmult}
\alias{\%*\%.rv}
\alias{matmult}
\title{Random Matrix Multiplication}
\description{
  Multiplies two random matrices, if they are conformable.
  If one argument is a vector, it will be coerced to either a row or
  column matrix to make the two arguments conformable.
  If both are vectors it will return the inner product.
}
\usage{
  \method{\%*\%}{rv}(x, y)
}
\arguments{
  \item{x, y}{numeric or complex matrices or vectors.}
}
\details{
  Optimized internally for the case of random matrix 
  multiplied by a constant column vector.
}
\value{
  The (distribution of the) matrix product.
  Use \code{\link{drop}} to get rid of dimensions
  which have only one level.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\seealso{
  \code{\link{matrix}},
  \code{\link{Ops}},
  \code{\link{diag}}.
}
\examples{
x <- 1:4
(z <- x \%*\% x)    # scalar ("inner") product (1 x 1 matrix)
drop(z)             # as scalar

y <- diag(x)
z <- matrix(1:12, ncol = 3, nrow = 4)
y \%*\% z
y \%*\% x
x \%*\% z
}
\keyword{array}
\keyword{arith}
