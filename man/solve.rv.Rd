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
Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations 
Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{classes}

