\name{rvdirichlet}
\title{Generate Random Variables from a Dirichlet Sampling Model}
\alias{rvdirichlet}
\usage{
  rvdirichlet(n=1, alpha)
}
\description{
  Generates random variables from a Dirichlet sampling model.
}
\arguments{
  \item{n}{integer: number of vectors to generate}
  \item{alpha}{the parameter vector; may be random}
}
\details{
  The Dirichlet distribution is a generalization of the Beta distribution.
  (If alpha is of length two, \code{rvdirichlet} draws from the Beta model.)

}
\value{
  A random vector (rv object) of length \code{n}.
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
  a <- rvdirichlet(1, alpha=c(6, 3, 1)) # 
  sum(a) # one with probability 1
}

\keyword{classes}

