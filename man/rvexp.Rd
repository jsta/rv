\name{rvexp}
\title{Generate Random Vectors from an Exponential Sampling Model}
\alias{rvexp}
\usage{
  rvexp(n=1, rate = 1)
}
\description{
  \code{rvexp}
}
\arguments{
  \item{n}{integer: number of variables to generate}
  \item{rate}{prior distribution for the rate parameter (constant or random)}
}
\details{
  \code{rvexp} 
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
  y <- rvexp(1, rate=rvexp(1)) # What marginal distribution does y have now?
}

\keyword{classes}

