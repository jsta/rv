\name{rvattr}
\title{Attributes of Random Variables}
\alias{rvattr}
\alias{rvattr<-}
\usage{
rvattr(x, attrib=NULL)
rvattr(x, attrib=NULL, by.name=FALSE) <- value
}
\description{
  \code{rvattr}
}
\arguments{
  \item{x}{an object}
  \item{attrib}{name of the attribute}
  \item{by.name}{logical; attempt matching of attributes by name?}
  \item{value}{vector of values to set; can be a list or an atomic vector}
}
\details{
  If \code{by.name=TRUE}, the values within the list \code{value} 
  are matched by their name (e.g. `\code{theta[1]}') if possible.
  Matching by \code{NA} or the empty string in a name is not possible.
  
  Otherwise, the list is matched by position; in this case, the length
  of \code{value} must be equal to that of \code{x}.
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
  ##
}

\keyword{classes}

