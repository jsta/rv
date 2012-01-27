\name{apply.rv}
\title{Apply Functions over Margins of Random Arrays}
\alias{apply.rv}
\usage{
  apply.rv(X, MARGIN, FUN, \dots)
}
\description{
  The \code{rv}-compatible version of \code{apply} 
}
\arguments{
 \item{X}{a random array}
 \item{MARGIN}{subscripts.}
 \item{FUN}{function.}
  \item{\dots}{optional arguments to \code{FUN}.}
}
\details{
  This is the rv-compatible version of the function \code{\link{apply}}.
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
\seealso{
  \code{\link{apply}}
}
\examples{
\dontrun{
  x <- rvmatrix(rvnorm(12), nrow=3, ncol=4)
  print(apply.rv(x, 1, sum))
}
}

\keyword{manip}

