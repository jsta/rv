\name{rvchisq}
\title{Generate Random Variables from a Chi-Square Sampling Model}
\alias{rvchisq}
\usage{
  rvchisq(n=1, df, ncp = 0)
}
\description{
  Generates a random vector from a chi-square sampling model.
}
\arguments{
 \item{n}{number of variables to generate}
 \item{df}{integer, degrees of freedom, may be random}
 \item{ncp}{non-centrality parameter, may be random}
}
\details{
If any of the arguments are random, 
the resulting simulations may have non-Poisson marginal distributions.
}
\value{
  A random vector (rv object) of length \code{n}.
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

