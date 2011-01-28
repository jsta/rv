\name{rvgamma}
\title{Generate Random Variables from a Gamma Sampling Model}
\alias{rvgamma}
\alias{rvngamma}
\usage{
  rvgamma(n=1, shape, rate = 1, scale = 1/rate)
  rvngamma(n=1, shape, rate = 1, scale = 1/rate)
}
\description{
  Generates random variables from a Gamma sampling model.
}
\arguments{
  \item{n}{integer: number of variables to generate}
  \item{shape}{shape parameter, may be a rv}
  \item{rate}{rate parameter, may be a rv}
  \item{scale}{inverse of rate, may be specified optionally instead of rate}
}
\details{
  \code{rvngamma(n, shape, rate)} is equivalent to
  \code{rvgamma(n, 1/3 + shape, rate)}.  
}
\value{
  A random vector (rv object).
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
  round(rvmedian(rvngamma(n=1, shape=1:10, rate=1)), 1) ## close to 1:10
}

\keyword{classes}

