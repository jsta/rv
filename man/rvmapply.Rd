\name{rvmapply}
\title{Apply a function to multiple random vector arguments}
\alias{rvmapply}
\alias{rvVectorize}
\usage{
  rvmapply(FUN, \dots, MoreArgs=NULL, SIMPLIFY = FALSE, USE.NAMES=TRUE, SAMPLESIZE=NULL)
  rvVectorize(FUN, vectorize.args = arg.names, SIMPLIFY = FALSE, USE.NAMES = TRUE, SAMPLESIZE=NULL) 
}
\description{
  \code{rvmapply} is the rv-compatible version of \code{\link{mapply}}.
  It repeats the function \code{FUN} for each joint draw of the 
  random (or constant) arguments, while allowing vectorizing.
}
\arguments{
  \item{FUN}{the function to apply to the simulations of \code{X}.}
  \item{MoreArgs}{Other args passed to \code{FUN} `as is' (must not be rv objects unless the function already accepts them)}
  \item{USE.NAMES}{logical; see \code{\link{mapply}} for details}
  \item{SIMPLIFY}{logical; see \code{\link{mapply}} for details}
  \item{SAMPLESIZE}{if specified, takes a (joint) sample of the simulations and processes only them.}
  \item{vectorize.args}{a character vector of arguments which should be vectorized. Defaults to all arguments to FUN.}
  \item{\dots}{further arguments to \code{FUN}, possibly random vectors or array.}
}
\details{
  \code{rvmapply} applies a given function to each 
  simulation (vector or array) of the given random vectors,
  returning a the results as a random vector or array. 

  The dimensions of each joint draw are preserved.
  For an example, see \code{\link{solve}}, that returns
  the distribution of the inverse of a random matrix.

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
}

\keyword{manip}

