\name{rvbeta}
\title{Generate Random Vectors from a Beta Sampling Model}
\alias{rvbeta}
\usage{
  rvbeta(n=1, shape1, shape2)
}
\description{
  \code{rvbeta}
}
\arguments{
  \item{n}{integer, number of random variables to generate}
  \item{shape1}{positive number or rv, 1st shape parameter}
  \item{shape2}{positive number or rv, 2nd shape parameter}
}
\details{
  \code{rvbeta} 
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
   n <- 12         # sample size
   y <- (0:(n-1))  # observations
   a <- b <- 1     # a uniform prior
   rvbeta(1, shape1=a+y, shape2=b+n-y)
}

\keyword{classes}

