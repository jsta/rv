\name{simapply}
\title{Apply a Function to Rows of Random Vectors}
\alias{simapply}
\usage{
simapply(x, FUN, \dots)
}
\description{
  \code{simapply} applies a given function \code{FUN} to each row of the simulation matrix,
  returning an rv object.
}
\arguments{
  \item{x}{a rv object}
\item{FUN}{a function}
\item{\dots}{further arguments passed to \code{FUN}}
}
\details{
  \code{simapply} ...
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
  x <- rvnorm(10)
  simapply(x, mean) # Same result as that of mean(x).
}

\keyword{classes}

