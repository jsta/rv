\name{outer.rv}
\title{Outer Product of Random Arrays}
\alias{outer.rv}
\usage{
outer.rv(X, Y=NULL, FUN="*", \dots)
}
\description{
  \code{outer.rv}
}
\arguments{
  \item{X}{First argument for function \code{FUN}}
  \item{Y}{Second argument for function \code{FUN}; if missing, \code{X} is used instead}
  \item{FUN}{a function to use on the outer products; a character string or a function}
  \item{\dots}{optional arguments to be passed to \code{FUN}}
}
\details{
  Implements the outer product for random arrays.

  See the function \code{outer} for further details.
}
\value{
  A random array.
}
\references{ 
Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations 
Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{classes}

