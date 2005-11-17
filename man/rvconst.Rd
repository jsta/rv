\name{rvconst}
\title{Random Vector with a Point-Mass Distribution}
\alias{rvconst}
\usage{
rvconst(n=1,x)
}
\description{
  Coerces a given vector of constants into a random vector
  with 1 simulation in each component.
}
\arguments{
  \item{n}{integer: number of variables to generate}
  \item{x}{a vector of constants}
}
\details{
  Coerces a given vector of constants into a random vector
  with 1 simulation in each component.
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

