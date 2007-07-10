\name{rep.rv}
\title{Replicate Elements of Random Vectors}
\alias{rep.rv}
\usage{
  \method{rep}{rv}(x, times, \dots)
}
\description{
  Transpose a random array by permuting its dimensions and optionally resizing it.
}
\arguments{
 \item{x}{a random vector to be replicated}
 \item{times}{number of replications}
 \item{\dots}{further arguments passed to \code{rep}}
}
\details{
 This is the rv-compatible version of the function \code{\link{rep}}.

 Since \code{rep} is not a generic function, 
 the whole name \code{rep.rv} must be specified when calling the function
 when \code{x} is an 'rv' object.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\seealso{
  \code{\link{rep}}
}
\examples{
  print(rep(rvnorm(1), times=4))
}

\keyword{manip}

