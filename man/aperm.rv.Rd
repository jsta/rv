\name{aperm.rv}
\title{Random Array Transposition}
\alias{aperm.rv}
\usage{
  aperm.rv(x, \dots)
}
\description{
  Transpose a random array by permuting its dimensions and optionally resizing it.
}
\arguments{
 \item{x}{the random matrix to be transposed}
  \item{\dots}{further arguments passed to \code{aperm}}
}
\details{
  This is the rv-compatible version of the function \code{\link{aperm}}.

  Since \code{aperm} is not a generic function, 
  the whole name \code{aperm.rv} must be specified when calling the function
  when \code{X} is an 'rv' object.
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
  \code{\link{aperm}}
}
\examples{
  x <- rvarray(rvnorm(24), dim=c(2,3,4))
  print(aperm.rv(x))
}

\keyword{manip}

