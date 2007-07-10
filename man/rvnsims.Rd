\name{rvnsims}
\title{Number of simulations stored in each component of an rv object}
\alias{rvnsims}
\usage{
  rvnsims(x)
}
\description{
  \code{rvnsims} returns the number of simulations 
  stored in each component of its argument.
}
\arguments{
  \item{x}{an rv object}
}
\details{
  If the argument is a number, \code{rvnsims} returns 1.
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
  rvnsims(1.23)             # one
  x <- rvnorm(1)            # an rv
  rvnsims(x)                # equal to setnsims()
  rvnsims(x)==nrow(sims(x)) # TRUE
  rvnsims(x)==setnsims()    # TRUE
}

\keyword{classes}

