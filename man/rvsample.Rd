\name{rvsample}
\title{Sample from the Simulation Matrix of a Random Variable}
\alias{rvsample}
\usage{
  rvsample(x, size=1, replace=TRUE, prob=NULL)
}
\description{
  Draws a sample of desired size from each component of a given random variable \code{x}.
}
\arguments{
  \item{x}{an object}
  \item{size}{size of the sample}
  \item{replace}{sample with replacement?}
  \item{prob}{optional probability vector, usually not needed}
}
\details{
  
}
\value{
  A \emph{numeric} array of dimensions \code{size} times \code{length(x)}.
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

