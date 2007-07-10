\name{rvsample}
\title{Draw a Sample from the Simulation Matrix of a Random Variable}
\alias{rvsample}
\usage{
  rvsample (x, size = 1, jointly = TRUE, reject.na = FALSE) 
}
\description{
  Draws a sample of desired size from each component of a given random variable \code{x}.
}
\arguments{
  \item{x}{an object}
  \item{size}{size of the sample}
  \item{jointly}{return joint simulations and not simulations from each component separately}
  \item{reject.na}{reject each draw that contains an \code{NA}}
}
\details{
  Samples (with replacement) from the distribution of the random variable object.
  In effect it samples from the rows of the simulation matrix \code{sims(x)}.
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

