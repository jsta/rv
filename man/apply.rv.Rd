\name{apply.rv}
\title{Distribution of the Range of a Random Vector}
\alias{apply.rv}
\usage{
  \method{apply}{rv}(X, MARGIN, FUN, \dots)
}
\description{
  \code{apply.rv}
  returns a 2-component random vector containing 
  the distributions of the minimum and the maximum values of all the given arguments.
}
\arguments{
 \item{X}{}
 \item{MARGIN}{}
 \item{FUN}{}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  This is the rv-compatible version of the function \code{\link{apply}}.
  Since \code{apply} is not a generic function, 
  the whole name \code{apply.rv} must be specified when calling the function
  when \code{X} is an 'rv' object.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\seealso{
  \code{\link{quantile.rv}}
}
\examples{
  x <- rvmatrix(rvnorm(12), nrow=3, ncol=4)
  print(apply.rv(x, 1, sum))
}

\keyword{manip}

