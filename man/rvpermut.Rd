\name{rvpermut}
\title{Random Vectors with a Permutation Distribution}
\alias{rvpermut}
\alias{rvpermut}
\usage{
  rvpermut(data, prob=NULL)
}
\description{
  Generates a random vector with each component having a permutation distribution
  based on the given (fixed) data vector.
}
\arguments{
  \item{data}{a fixed numeric vector}
  \item{prob}{optional probabilities for the components in \code{data}}
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
  x <- rvpermut(1:10)
}

\keyword{classes}

