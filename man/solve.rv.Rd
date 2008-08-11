\name{solve.rv}
\title{Random Vectors}
\alias{solve.rv}
\usage{
  \method{solve}{rv}(a, b, \dots)
}
\description{
  \code{solve.rv}
}
\arguments{
 \item{a}{a square random vector containing the coefficients of the linear system}
 \item{b}{a square random vector giving the right-hand side(s) of the linear system}
 \item{\dots}{further arguments passed to \code{solve}}
}
\details{
  \code{solve.rv} is the rv-object compatible version of the function \code{solve}.

  For details of the function, see \code{\link{solve}}.
}
\seealso{
  \code{\link{solve}} 
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
  #
}

\keyword{classes}

