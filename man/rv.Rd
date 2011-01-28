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
\alias{as.rvobj}
\alias{is.rvobj}
\usage{
rv(length = 0)
as.rv(x, \dots)
is.rv(x)
is.random(x)
as.rvobj(x)
is.rvobj(x)
}
\description{
  Creates or tests for objects of type ``\code{rv}''.
}
\arguments{
  \item{length}{desired length.}
  \item{x}{object to be coerced or tested.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  \code{rv} creates a random vector of the specified length. 
The elements of the vector are all equal to \code{NA}.

  \code{is.rv} returns TRUE if its argument is a rv object, FALSE otherwise.

  \code{as.rv} attempts to coerce its argument to the random vector (rv) type.

  \code{is.random} returns \code{TRUE} or \code{FALSE}
  for each component of the argument vector, depending on whether the component is 
  a random variable object.

  \code{is.rvobj} tests whethe its argument object
  is either of class \code{rv} or of class \code{rvsummary}.

  \code{as.rvobj} coerces its argument object to \code{rv} unless
  the object is an rv object (\code{is.rvobj(x)} is \code{TRUE}).
}
\value{
  An rv object of desired length, with the single simulation value \code{NA}.

}
\note{
  rv objects are internally lists with the class attribute set to ``\code{rv}".
  The number of simulations in rv objects is set by \code{\link{setnsims}}.
  This is by default set to 2500.
}
\seealso{
  For a short version of the paper, view the vignette by \code{vignette("rv")}.
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
  x <- rv(1)
  
}
\keyword{classes}

