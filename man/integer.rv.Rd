\name{as.integer.rv}
\title{Integer Random vectors}
\alias{as.integer.rv}
\usage{
  \method{as.integer}{rv}(x, \dots)
}
\description{
  Coerces a random variable to an integer-valued (discrete) one
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{Further arguments passed on}
}
\details{
   In effect, the function  \code{as.integer} is applied to all simulations.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\seealso{
  \code{\link{as.logical.rv}}.
}
\note{
  \code{is.integer(x)} returns \code{TRUE} if and only if \emph{each} component 
  of \code{x} is integer-valued (each simulation vector is of type 'integer').
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  x <- rvpois(lambda=3)   # some integer-valued random variable
  print(x)
  is.integer(x)           # FALSE, because by default x is 'double'!
  x <- as.integer(x)      # coerce to integer
  is.integer(x)           # TRUE
  print(x)                # Shows also the 'min' and 'max' columns
}

\keyword{classes}

