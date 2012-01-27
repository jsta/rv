\name{Covariance, Correlation, Variance, and Standard deviation Statistics}
\title{Distributions of Various Statistics of a Random Vector and Array}
\alias{cov.rv}
\alias{cor.rv}
\alias{var.rv}
\alias{sd.rv}
\usage{
  sd.rv(x, na.rm=FALSE)
  cov.rv(x, y=NULL, \dots)
  cor.rv(x, y=NULL, \dots)
  var.rv(x, \dots)
}
\description{
 \code{var.rv}, \code{cov.rv} and \code{cor.rv} compute the distribution
 of the variance statistic of x and the distribution of the covariance
 statistic or the correlation statistic of x and y if these are vectors.
 If x and y are matrices then the covariances (or correlations)
 between the columns of x and the columns of y are computed.
}
\arguments{
  \item{x}{a numeric or random vector, matrix, or a data frame}
  \item{y}{\code{NULL} (default) or a vector, matrix or data frame with compatible dimensions to x. The default is equivalent to y = x (but more efficient).}
  \item{\dots}{further arguments passed to the corresponding numeric functions}
}
\details{

  These functions are compatible with \emph{both} numeric and rv
  objects.
  To make your code compatible with \code{rv} objects, use
  e.g. \code{sd.rv} instead of \code{sd}. 

  The functions \code{cov.rv} is implemented by
  applying the corresponding numerical  
  function to the rows of the simulation matrices
  of \code{x} and \code{y} and forming a new \code{rv} object
  from the resulting vector of simulations.
  Alternatively \code{x} may be a random matrix (and \code{y} \code{NULL}).
  %Then the numerical function \code{cov}.
  
  \code{cor.rv} works similarly, but returns the distribution of the
  correlation statistic (i.e. function).

  \code{var.rv} computes the distribution of the variance statistic.
  \code{sd.rv} is the square root of the result obtained by \code{var.rv}.
}
\value{
  A random vector or array.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
%\seealso{
  %\code{\link[pkg:stats]{cov}}, \code{\link[pkg:stats]{cor}}, \code{\link[pkg:stats]{var}}, \code{\link[pkg:stats]{sd}} 
  %for details.
%}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  #
}

\keyword{internal}

