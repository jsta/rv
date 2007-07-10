\name{atomic}
\title{Random Vectors Are Atomic}
\alias{atomic}
\alias{is.atomic}
\usage{
  is.atomic(x)
}
\description{
  Tests whether an object is ``atomic".
  For compatibility, 'rv' objects are defined to be atomic. 
}
\arguments{
  \item{x}{an object, random or constant}
}
\details{
  If \code{x} is not an rv object, the function returns the value of
  the primitive function \code{"is.atomic"}.

  The original function \code{is.atomic} is replaced with
  the function in the rv package; the original function is 
  restored upon detaching the package.
}
\value{
  \code{TRUE} or \code{FALSE}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  is.atomic(rvnorm(1))
}

\keyword{classes}

