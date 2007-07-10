\name{rvsims}
\title{Create a Random Vector from a Matrix of Simulations}
\alias{rvsims}
\usage{
  rvsims(sims, n.sims=setnsims(), permute=FALSE, save.order=FALSE)
}
\description{
\code{rvsims} takes a vector, matrix, or 3-way array (\code{sims})
containing simulations,
and returns a random vector (an object of type 'rv')
}
\arguments{
\item{sims}{an array of simulations (1, 2, or 3-dimensional)}
\item{n.sims}{number of simulations to save}
\item{permute}{logical, indicate if scramble the simulations}
\item{save.order}{logical, indicate if the original order of simulations should be saved (only if permute is \code{TRUE})}
}
\details{
A vector is interpreted to contain simulations for one single random variable;
if \code{sims} is two-dimensional, the columns are supposed to
contain simulations for several random variables.
If \code{sims} is three-dimensional, it is supposed to be output from
a Markov chain simulation program:
the first dimension corresponds to the simulations, 
the second corresponds to the chains (typically 3 or 5),
and the last dimension corresponds to the variables.

If \code{permute} is \code{TRUE}, the simulations are scrambled,
i.e. the \emph{rows} are permuted randomly.
If \code{permute=TRUE} \emph{and} \code{save.order=TRUE},
the original order of the simulations
is saved in an attribute \code{rvsim$order}.

If \code{sims} is three-dimensional, the simulations are always scrambled.

}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  # x and y have the same distributions:
  x <- sims(rnorm(1000))
  y <- rvnorm(1)
}

\keyword{classes}

