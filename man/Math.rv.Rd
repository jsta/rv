\name{Math.rv}
\title{Mathematical functions and Operators for rv Objects}
\alias{Math.rv}
\alias{Ops.rv}
\alias{!.rv}
\alias{Math.rvsim}
\alias{Ops.rvsim}
\alias{cumsum.rv}
\alias{cumprod.rv}
\alias{cummin.rv}
\alias{cummax.rv}
\usage{
  Math.rv(x, ...)
  Ops.rv(e1, e2)
}
\description{
  Mathematical functions and operators adapted to work with 
  random variable (rv) objects.
}
\arguments{
  \item{x}{object}
  \item{e1}{object}
  \item{e2}{object}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvnorm(10)
  -x
  x + 1:10
  cumsum(x)
  cumprod(exp(x))
}

\keyword{classes}
