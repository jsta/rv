\name{rvt}
\title{Generate Random Variables from a Student-t Sampling Model}
\alias{rvt}
\usage{
  rvt(n=1, mu=0, scale=1, df, ncp, Sigma)
}
\description{
  Generates a random variable from a Student-t sampling model.
}
\arguments{
  \item{n}{integer, number of scalars to generate}
  \item{mu}{location, may be a rv}
  \item{scale}{scale, may be a rv}
  \item{ncp}{non-centrality parameter}
  \item{df}{degrees of freedom, may be a rv}
  \item{Sigma}{(optional) scaling matrix for multivariate generation}
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
  df <- 3
  x <- rvt(n=1, df=df)
  y <- rvnorm(1)/sqrt(rvchisq(1, df=df)/df) # Same distribution as above
  print(c(x,y))
}

\keyword{classes}

