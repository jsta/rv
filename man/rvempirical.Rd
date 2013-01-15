\name{rvempirical}
\title{Generate a Random Vector from an Empirical Distribution}
\alias{rvempirical}
\usage{
  rvempirical(n, data)
}
\description{
  \code{rvempirical} generates a random vector of the same length as data
  from the empirical distribution of the data.
}
\arguments{
  \item{n}{Number of i.i.d. rv components to generate}
  \item{data}{Data (constants)}
}
\details{
  \code{rvempirical} 
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
  y <- c(1.0, 1.2, 3, 1.1, 0.8, 0.9) ## Some data
  x <- rvempirical(4, data=y) # A random vector of length 4 (each component has the same distribution)
}

\keyword{classes}

