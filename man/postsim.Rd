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
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- 1:20
  y <- rnorm(length(x), mean=x, sd=10)
  print(summary(lm.fit <- lm(y ~ x)))
  random.estimates <- postsim(lm.fit)
}

\keyword{classes}


