\name{Summary.rv}
\title{Summary}
\alias{Summary.rv}
\usage{
  \method{Summary}{rv}(\dots, na.rm=FALSE)
}
\description{
  \code{Summary.rv}
}
\arguments{
  \item{\dots}{rv object}
  \item{na.rm}{logical, remove NAs?}
}
\details{
    
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
  x <- rvnorm(5)
  all(x>0)
  any(x>2)
  sum(x)
  min(x)
  max(x)
}

\keyword{internal}

