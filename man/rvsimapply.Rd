\name{rvsimapply}
\title{Apply a Function to Columns of the Matrix of Simulation of a Random Vector}
\alias{rvsimapply}
\usage{
rvsimapply(x, FUN, \dots)
}
\description{
  \code{rvsimapply}
}
\arguments{
  \item{x}{an object}
  \item{FUN}{}
  \item{\dots}{further arguments passed to the function \code{FUN}}
}
\details{
  \code{rvsimapply} applies a given function to 
  the \emph{rows} of the simulation matrix of the given random vector.

  If the function is to be applied to \emph{rows} of the simulation matrix,
  use \code{\link{simapply}} or \code{\link{mapply.rv}} instead.

  Usually used in functions that implement an 'rv'-compatible routine.  
}
\value{
  A numeric vector or array.
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

