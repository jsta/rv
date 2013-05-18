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
If \code{sims} is a plain numeric vector, this is interpreted to be
equivalent to a one-dimensional array, containing simulations for one
single random variable.  

If the array \code{sims} is one-dimensional, this is interpreted to
be equivalent to a two-dimensional array with 1 column.

If \code{sims} is two-dimensional, the \emph{columns} are supposed to
contain simulations for one or more several random variables.

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
  ## x and y have the same distributions but not the same simulations:
  n.sims <- 200L
  setnsims(n.sims)
  y <- rvnorm(1)
  x1 <- rvsims(rnorm(n.sims))
  ##
  s <- sims(x1)
  z <- array(s) ## One-dimensional array
  x2 <- rvsims(z) ## Same as 
  ##
  identical(x1, x2) ## TRUE
  ##
  s <- t(array(rnorm(n.sims * 2, mean=c(0, 10)), dim=c(2, n.sims)))
  x3 <- rvsims(s)
  identical(2L, length(x3)) ## TRUE
}

\keyword{classes}

