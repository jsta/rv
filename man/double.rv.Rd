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
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
%\seealso{
%  
%}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  x <- as.logical(rvbern(prob=0.5)) 
  print(x)
  print(as.double(x))
}

\keyword{classes}

