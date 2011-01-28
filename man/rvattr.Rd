\name{rvattr}
\title{Attributes of Random Variables}
\alias{rvattr}
\alias{rvattr<-}
\usage{
rvattr(x, attrib=NULL)
rvattr(x, attrib=NULL) <- value
}
\description{
  \code{rvattr}
}
\arguments{
  \item{x}{an object}
  \item{attrib}{name of the attribute}
  \item{value}{list of values to set}
}
\details{
  The values within the list \code{value} 
  are matched by their name (e.g. `\code{theta[1]}') if possible;
  if both \code{x} and \code{value} do not have names,
  the list is matched by position. 
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

