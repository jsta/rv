\name{rvexp}
\title{Generate Random Vectors from an Exponential Sampling Model}
\alias{rvexp}
\usage{
  rvexp(n=1, rate = 1)
}
\description{
  \code{rvexp}
}
\arguments{
  \item{n}{integer: number of variables to generate}
  \item{rate}{prior distribution for the rate parameter (constant or random)}
}
\details{
  \code{rvexp} 
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  y <- rvexp(1, rate=rvexp(1)) # What marginal distribution does y have now?
}

\keyword{classes}

