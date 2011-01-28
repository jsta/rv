\name{aperm.rv}
\title{Random Array Transposition}
\alias{aperm.rv}
\usage{
  \method{aperm}{rv}(a, perm, \dots)
}
\description{
  Transpose a random array by permuting its dimensions and optionally resizing it.
}
\arguments{
  \item{a}{the random matrix to be transposed}
  \item{perm}{the subscript permutation vector. See the manual page for
    the gneric method aperm.}
  \item{\dots}{further arguments passed to \code{aperm}}
}
\details{
  This is the rv-compatible version of the function \code{\link{aperm}}.
  It first applies 
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
  print(aperm(x))
}

\keyword{manip}

