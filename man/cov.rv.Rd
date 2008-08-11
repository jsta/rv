\name{Covariance, Correlation, Variance, and Standard deviation Statistics}
\title{Distributions of Various Statistics of a Random Vector and Array}
\alias{cov}
\alias{cor}
\alias{var}
\alias{sd}
\usage{
  cov(x, \dots)
  cor(x, \dots)
  var(x, \dots)

  \method{cov}{rv}(x, y=NULL, \dots)
  \method{cor}{rv}(x, y=NULL, \dots)
  \method{var}{rv}(x, \dots)

  sd(x, na.rm=FALSE)

}
\description{
 \code{var}, \code{cov} and \code{cor} compute the distribution of the variance statistic of x and the distribution of the covariance statistic or the correlation statistic of x and y if these are vectors. If x and y are matrices then the covariances (or correlations) between the columns of x and the columns of y are computed.
}
\arguments{
  \item{x}{a numeric or random vector, matrix, or a data frame}
  \item{y}{\code{NULL} (default) or a vector, matrix or data frame with compatible dimensions to x. The default is equivalent to y = x (but more efficient).}
  \item{\dots}{further arguments passed to corresponding numeric functions}
}
\details{

  The functions \code{cov}, \code{var}, \code{cor}
  are now generic - whenever the \code{rv} package is loaded.

  The original functions in the \code{stats} package 
  are the default functions of these generic functions.

  To access the original help page in the stats package,
  go to \code{\link[pkg:stats]{cov2cor}} or
  detach the rv package and then access the help pages in the stat package.

  \code{sd} is not generic but it is rewritten to catch the possibility 
  that the argument is a random variable object.

  \code{cov} gives the distribution (that is, a random variable object)
  of the covariance statistic, as computed by the 
  standard R function \code{\link[pkg:stats]{cov}}. 

  This is implemented simply by applying the numerical 
  \code{cov} function to the rows of the simulation matrices
  of \code{x} and \code{y} and forming a new rv object
  from the resulting vector of simulations.
  Alternatively \code{x} may be a random matrix (and \code{y} \code{NULL}).
  Then the numerical function \code{cov} is m
  
  \code{cor} works similarly, but returns the distribution of the
  correlation statistic (i.e. function).

  \code{var} computes the distribution of the variance statistic.

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
\seealso{
  \code{\link[pkg:stats]{cov}}, \code{\link[pkg:stats]{cor}}, \code{\link[pkg:stats]{var}}, \code{\link[pkg:stats]{sd}} 
  for details.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  #
}

\keyword{internal}

