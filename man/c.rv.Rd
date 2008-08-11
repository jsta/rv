\name{c}
\title{Concatenation of Random Vectors}
\alias{c.rv}
\usage{
  \method{c}{rv}(\dots, recursive = FALSE)
}
\description{
  Concatentates random vectors.
}
\arguments{
  \item{\dots}{objects to be concatenated. Can be a mixture of constants and rv objects.}
  \item{recursive}{logical. If recursive = TRUE, the function recursively descends through lists (and pairlists) combining all their elements into a vector.}
}
\details{
  NOTE: \code{recursive} has not yet been tested.    
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
  x <- rvnorm(2)
  y <- rvbern(2, prob=0.5)
  z <- c(x, y)
  print(z)
  z1 <- c(1, z)
  z2 <- c(as.rv(1), z)
  print(z1)
  print(z2)
}

\keyword{classes}

