\name{sort}
\title{Distribution of Order Statistics of a Random Vector}
\alias{sort}
\usage{
  sort(x, \dots)
}
\description{
  \code{sort} is a replacement for the standard \code{\link[base]{sort}}
  function, which is called if \code{x} is not an \code{rv} object.
}
\arguments{
  \item{x}{an object}
  \item{\dots}{further arguments passed to \code{sort}}
}
\details{
  Once \code{sort} becomes a generic (and truly object-oriented) method,
  this function will also become a truly object-oriented method.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\seealso{
  The original help page for  \code{\link[base]{sort}}.
  Other adapted functions: \link{rvAdapted}, \code{\link{cov}}.
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

