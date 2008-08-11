\name{constant}
\title{Constant Vectors}
\alias{is.constant}
\alias{as.constant}
\alias{as.constant.rv}
\usage{
 is.constant(x)
 as.constant(x)
 \method{as.constant}{rv}(x)
}
\description{
  Tests or coerces objects that are non-random.  
}
\arguments{
  \item{x}{an object, random variable (rv) or not}
}
\details{
  \code{is.constant} returns \code{TRUE}
  for each component of the argument object
  if there is only one simulation (that is, the variable is ``constant").

  Note: rv objects that merely have variance zero are not therefore necessarily 
  constants. 
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
  is.constant(1)         # TRUE
  is.constant(as.rv(1))  # TRUE
  setnsims(200)
  x <- rvbern(prob=0.001)
  all(sims(x)==0)        # most probably true
  is.constant(x)         # always FALSE
}

\keyword{classes}

