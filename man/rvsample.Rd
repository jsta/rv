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
  #
}

\keyword{classes}

