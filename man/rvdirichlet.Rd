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
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  a <- rvdirichlet(1, alpha=c(6, 3, 1)) # 
  sum(a) # one with probability 1
}

\keyword{classes}

