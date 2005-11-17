\name{constant}
\title{Constant Vectors}
\alias{is.constant}
\alias{is.constant.rv}
\alias{is.constant.default}
\alias{is.constant.rvsim}
\alias{as.constant}
\alias{as.constant.rvsim}
\alias{as.constant.rv}
\usage{
 is.constant(x)
 \method{is.constant}{rv}(x)
 \method{is.constant}{default}(x)
 \method{is.constant}{rvsim}(x)
 as.constant(x)
 \method{as.constant}{rvsim}(x)
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
  if the implied variance (sample variance of the simulations) is zero
  or if there is only one simulation (that is, the variable is ``constant").
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  is.constant(1)         # TRUE
  is.constant(rvnorm(1)) # FALSE
}

\keyword{classes}

%%% add as.constant.