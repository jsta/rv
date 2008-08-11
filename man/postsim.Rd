\name{postsim}
\title{Generate Posterior Simulations for lm or glm Objects}
\alias{postsim}
\alias{postsim.lm}
\alias{postsim.glm}
\usage{
  postsim(fit)
  \method{postsim}{lm}(fit)
  \method{postsim}{glm}(fit)
}
\description{
  Generate posterior simulations for a given fitted linear or general linear model,
  assuming the standard "noninformative" priors on the unknowns.
}
\arguments{
  \item{fit}{an lm or glm object}
}
\value{
  A (named) random vector for each fitted coefficient.
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
  bayes.estimates <- postsim(fit)
  }
}

\keyword{classes}


