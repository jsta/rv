\name{rv}
\title{Random Vectors}
\alias{rv}
\alias{is.rv}
\alias{is.rv.rv}
\alias{is.rv.default}
\alias{as.rv}
\alias{as.rv.rv}
\alias{as.rv.numeric}
\alias{as.rv.integer}
\alias{as.rv.logical}
\alias{as.rv.list}
\alias{as.rv.matrix}
\alias{as.rv.default}
\alias{is.random}
\usage{
rv(length = 0)
as.rv(x)
is.rv(x)
is.random(x)
}
\description{
  Creates or tests for objects of type ``\code{rv}".
}
\arguments{
  \item{length}{desired length.}
  \item{x}{object to be coerced or tested.}
%  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  \code{rv} creates a random vector vector of the specified length. The elements of the vector are all equal to 0.

  \code{is.rv} returns TRUE if its argument is a rv object, FALSE otherwise.

  \code{as.rv} attempts to coerce its argument to the random vector (rv) type.

  \code{is.random} returns \code{TRUE} or \code{FALSE}
  for each component of the argument vector, depending on whether the component is 
  a random variable object.
}
\note{
  rv objects are internally lists with the class attribute set to ``\code{rv}".
}
\seealso{
  For a short version of the paper, type \code{doc("rv")}.
}
\references{Kerman, Jouni and Gelman, Andrew. 
Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
\url{http://www.stat.columbia.edu/~kerman/Research/postsim.pdf}
}

\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}

\examples{
  x <- rv(1)
}
\keyword{classes}

