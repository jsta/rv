\name{constant}
\title{Constant Vectors}
\alias{is.constant}
\alias{as.constant}
\alias{as.constant.rv}
\alias{as.constant.rvsummary}
\usage{
 is.constant(x)
 as.constant(x)
 \method{as.constant}{rv}(x)
 \method{as.constant}{rvsummary}(x)
}
\description{
  Functions to coerce or test for non-random objects.
}
\arguments{
  \item{x}{an object, random variable (rv) or not}
}
\details{
  \code{is.constant} returns \code{TRUE}
  for each component of the argument object
  if there is only one simulation (that is, the variable is ``constant'').

  Note: rv objects that merely have variance zero are not therefore necessarily 
  ``true'' constants. 

  \code{as.constant} coerces \code{rv} or \code{rvsummary} objects into 
  constant strings; \code{NA} is returned for component that is not random.
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
  x <- rvnorm(3)
  x[1] <- 1
  as.constant(x)         # 1, NA, NA
  all(is.random(x) & is.na(as.constant(x))) # always TRUE
}

\keyword{classes}

