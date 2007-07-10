\name{as.double.rv}
\title{Coercing Random Vectors to Real-valued}
\alias{as.double.rv}
\alias{as.real.rv}
\usage{
  \method{as.double}{rv}(x, \dots)
  \method{as.real}{rv}(x, \dots)
}
\description{
  Coerces random vector objects into double-valued ones.
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{other arguments}
}
\details{
   \code{as.double} coerces an rv object into double-valued one.
   In effect, the function  \code{as.double} is applied to all simulations.

   \code{as.real} is a synonym.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\seealso{
  \code{\link{integer.rv}}
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- as.logical(rvbern(prob=0.5)) 
  print(x)
  print(as.double(x))
}

\keyword{classes}

