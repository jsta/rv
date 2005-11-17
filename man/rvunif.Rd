\name{rvunif}
\title{Generate Random Vectors from a Uniform Sampling Model}
\alias{rvunif}
\usage{
rvunif(n=1, min=0, max=1)
}
\description{
  Generates random variables from a Uniform sampling model.
}
\arguments{
  \item{n}{integer: number of scalars to generate}
  \item{min}{lower limit of the distribution, (may be random)}
  \item{max}{upper limit of the distribution, (may be random)}
}
\details{
  
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  y <- rvunif(1, min=rvunif(1)-1, rvunif(1)+1) # What marginal distribution does y have now?
}

\keyword{classes}

