\name{simapply}
\title{Apply a Function to Rows of Simulations of Random Vectors}
\alias{simapply}
\usage{
simapply(x, FUN, \dots)
}
\description{
  \code{simapply} applies a given function \code{FUN} to each row of the simulation matrix,
  returning an rv object.
}
\arguments{
  \item{x}{a random vector}
 \item{FUN}{a function}
 \item{\dots}{further arguments passed to \code{FUN}}
}
\details{
  \code{simapply} applies a given function to 
  the \emph{rows} of the simulation matrix of the given random vector.

  If the function accepts \emph{arrays}, use \code{\link{mapply.rv}} instead.

  If the function is to be applied to \emph{columns} of the simulation matrix,
  use \code{\link{rvsimapply}} instead.

  Usually used in functions that implement an 'rv'-compatible routine.  
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\value{
  
  
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
  x <- rvnorm(10)
  simapply(x, mean) # Same result as that of mean(x).
}

\keyword{manip}

