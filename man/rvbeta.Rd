\name{rvbeta}
\title{Generate Random Vectors from a Beta Sampling Model}
\alias{rvbeta}
\alias{rvnbeta}
\usage{
  rvbeta(n=1, shape1, shape2)
  rvnbeta(n=1, shape1, shape2)
}
\description{
  \code{rvbeta} generates a random vector from the beta sampling model;

  \code{rvnbeta(n, a, b)} ("neutral" Beta distribution) is equivalent to
  \code{rvbeta(n, 1/3+a, 1/3+b)}.
}
\arguments{
  \item{n}{integer, number of random variables to generate}
  \item{shape1}{positive number or rv, 1st shape parameter}
  \item{shape2}{positive number or rv, 2nd shape parameter}
}
%\details{
%  \code{rvbeta} 
%}
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
   a <- b <- 1/3   # the neutral beta prior
   rvbeta(1, shape1=a+y, shape2=b+n-y)
   rvnbeta(1, shape1=y, shape2=n-y)
}

\keyword{classes}

