\name{mapply.rv}
\title{Apply a function to multiple random vector arguments}
\alias{mapply.rv}
\usage{
  mapply.rv(FUN, \dots, MoreArgs=NULL, USE.NAMES=TRUE, SAMPLESIZE=NULL)
}
\description{
  \code{mapply.rv} is a rv-compatible version of \code{\link{mapply}}.
  It repeats the function \code{FUN} for each joint draw of the 
  random (or constant) arguments.
}
\arguments{
  \item{FUN}{the function to apply to the simulations of \code{X}.}
  \item{MoreArgs}{Other args passed to \code{FUN} `as is'}
  \item{USE.NAMES}{logical; see \code{\link{mapply}} for details}
  \item{SAMPLESIZE}{number of simulations to process. This should be not specified whenever the function is numeric.}
  \item{\dots}{further arguments to \code{FUN}, possibly random vectors or array.}
}
\details{
  \code{mapply.rv} applies a given function to each 
  simulation (vector or array) of the given random vectors,
  returning a the results as a random vector or array. 

  Usually used in functions that implement an 'rv'-compatible routine.

  For an example of a function that uses \code{SAMPLESIZE}, 
  \code{\link{abline}}.
}
\note{
  If the function (\code{FUN}) has an argument ``\code{FUN}",
  it must be specified within the list supplied to \code{MoreArgs}.
}
\value{
  Depends on \code{FUN}; a random vector or array if \code{FUN} is numeric.
}
\seealso{
  \code{\link{mapply}}, \code{\link{simapply}}
}
\references{ 
Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations 
Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{manip}

