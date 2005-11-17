\name{rvgamma}
\title{Generate Random Variables from a Gamma Sampling Model}
\alias{rvgamma}
\usage{
  rvgamma(n=1, shape, rate = 1, scale = 1/rate)
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
  \code{rvgamma} 
}
\value{
  A random vector (rv object).
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{classes}

