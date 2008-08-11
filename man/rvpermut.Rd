\name{rvpermut}
\title{Random Vectors with a Permutation Distribution}
\alias{rvpermut}
\alias{rvpermut}
\usage{
  rvpermut(data, prob=NULL)
}
\description{
  Generates a random vector with each component having a permutation distribution
  based on the given (fixed) data vector.
}
\arguments{
  \item{data}{a fixed numeric vector}
  \item{prob}{optional probabilities for the components in \code{data}}
}
\details{
        
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
  x <- rvpermut(1:10)
}

\keyword{classes}

