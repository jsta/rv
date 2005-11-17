\name{nsims}
\title{Number of simulations stored in each component of an rv object}
\alias{nsims}
\usage{
  nsims(x)
}
\description{
  \code{nsims} returns the number of simulations 
  stored in each component of its argument.
}
\arguments{
  \item{x}{an rv object}
}
\details{
  If the argument is a number, \code{nsims} returns 1.
}
\value{
  A vector of integers.
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
  nsims(1.23)
  x <- rvnorm(1)
  nsims(x)
  nsims(x)==nrow(sims(x)) # Should be TRUE
  nsims(x)==rvnsims()     # Should be TRUE
}

\keyword{classes}

