\name{rvsims}
\title{Create Random Vectors from Simulation Draws}
\alias{rvsims}
\usage{
  rvsims(sims, n.sims=getnsims(), permute=FALSE)
}
\description{
\code{rvsims} takes a vector, matrix, or list (\code{sims})
containing simulations,
and returns a random vector (an object of type 'rv')
}
\arguments{
\item{sims}{an array of simulations (1, or 2-dimensional) or a list}
\item{n.sims}{number of simulations to save}
\item{permute}{logical, indicate if scramble the simulations}
}
\details{
A vector is interpreted to contain simulations for one single random variable.

If \code{sims} is two-dimensional, the columns are supposed to
contain simulations for several random variables.

If \code{sims} is a list, the numeric vectors are recursively combined to
a list of random vectors: each component of the list is supposed to be
containing \emph{one} (joint) draw from some distribution---this may be a list.

If \code{permute} is \code{TRUE}, the simulations are scrambled,
i.e. the joint draws are permuted randomly.
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
  # x and y have the same distributions:
  setnsims(200)
  y <- rvnorm(1)
  x <- sims(rnorm(200))
  # 
}

\keyword{classes}

