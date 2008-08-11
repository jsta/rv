\name{unlist.rv}
\title{Flatten Lists Containing rv Objects}
\alias{unlist.rv}
\usage{
  \method{unlist}{rv}(x, recursive = TRUE, use.names = TRUE)
}
\description{
  Given a list structure \code{x}, \code{unlist} simplifies it to produce a vector which contains all the atomic components (\emph{containing rv objects}) which occur in \code{x}.
}
\arguments{
 \item{x}{An R object, typically a list or vector (containing rv objects)}
  \item{recursive}{logical. Should unlisting be applied to list components of x?}
  \item{use.names}{logical. Should names be preserved? (now fixed to TRUE)}
}
\details{
  This is the rv-compatible version of the function \code{\link{unlist}}.

  Since \code{unlist} is not a generic function, 
  the whole name \code{unlist.rv} must be specified when calling the function
  when \code{x} is an 'rv' object.
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
  \code{\link{unlist}}
}
\examples{
  x <- list(a=rvnorm(2), b=rvnorm(3))
  print(unlist.rv(x))
}

\keyword{manip}

