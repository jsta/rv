% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvpermut.R
\name{rvpermut}
\alias{rvpermut}
\title{Random Vectors with a Permutation Distribution}
\usage{
rvpermut(data, prob = NULL)
}
\arguments{
\item{data}{a fixed numeric vector}

\item{prob}{optional probabilities for the components in \code{data}}
}
\description{
Generates a random vector with each component having a permutation
distribution based on the given (fixed) data vector.
}
\examples{

  x <- rvpermut(1:10)

}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects. Statistics and
Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\author{
Jouni Kerman \email{jouni@kerman.com}
}
\keyword{classes}
