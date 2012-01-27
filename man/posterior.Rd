\name{posterior}
\title{Generate Posterior Simulations}
\alias{posterior}
\alias{posterior.lm}
\alias{posterior.glm}
\usage{
  posterior(obj, \dots)
  \method{posterior}{lm}(obj, \dots)
  \method{posterior}{glm}(obj, \dots)
}
\description{
  Generate posterior simulations for a given fitted linear or general linear model,
  assuming the standard "noninformative" priors on the unknowns.
}
\arguments{
  \item{obj}{an object}
  \item{\dots}{further arguments}
}
\value{
  A (named) list of random vectors.
  For example, the \code{lm} method returns
  a list with components \code{sigma} (the residual s.d.)
  and \code{beta}, the regression coefficients.
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
  \dontrun{
  x <- 1:20
  y <- rnorm(length(x), mean=x, sd=10)
  print(summary(fit <- lm(y ~ x)))
  bayes.estimates <- posterior(fit)
  }
}

\keyword{classes}


