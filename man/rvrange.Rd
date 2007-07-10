\name{rvrange}
\title{Range of Simulations in a Random Variable}
\alias{rvrange}
\usage{
  rvrange(x, \dots)
}
\description{
  Returns the range (min, max) of the simulations for each component of the 
  given random vector or array.
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{arguments passed to \code{range}}
}
\details{
  \code{rvrange} applies the function \code{range} 
  to each component of the 'rv' object.
}
\value{
  A numeric (constant) array with 2 rows (=min, max) and \code{length(x)} columns.
}
\seealso{
  \code{\link{range.rv}} returns the \emph{distribution} of the range of a random vector.
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

