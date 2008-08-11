\name{as.bugs.rv}
\title{Coerce an rv object to a bugs object}
\alias{as.bugs.rv}
\usage{
  as.bugs.rv (x, DIC = FALSE)
}
\description{
  Coerces an rv object to a bugs object.
}
\arguments{
  \item{x}{an rv}
  \item{DIC}{logical, if TRUE, compute deviance information criterion (DIC)}
}
\value{
  a bugs object.
}
\note{
  Works currently ONLY with random variable vectors that have
  been created from 3-dimensional arrays containing the MCMC chains.
  This may change in the future.
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
%\examples{
%
%}

\keyword{classes}

