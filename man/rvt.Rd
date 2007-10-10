\name{rvt}
\title{Generate Random Variables from a Student-t Sampling Model}
\alias{rvt}
\usage{
  rvt(n=1, mu=0, scale=1, df, Sigma)
}
\description{
  Generates a random variable from a Student-t sampling model.
}
\arguments{
  \item{n}{integer, number of scalars to generate}
  \item{mu}{location, may be a rv}
  \item{scale}{scale, may be a rv}
  \item{Sigma}{(optional) scaling matrix for multivariate generation}
  \item{df}{degrees of freedom, may be a rv}
}
\details{
  This function generates both univariate (independent and identically distributed)
  Student-t random variables and multivariate Student-t distributed vectors
  (with a given scaling matrix).

  For details of the parameters, see the entry on \code{mvt}
  in the \code{mvtnorm} package.
}
\note{
If any of the arguments are random, 
the resulting simulations may have non-t marginal distributions.
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

