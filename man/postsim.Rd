\name{postsim}
\title{Generate Posterior Simulations for lm or glm Objects (defunct)}
\alias{postsim}
\alias{postsim.lm}
\alias{postsim.glm}
\usage{
  postsim(fit)
}
\description{

  \emph{DEFUNCT.} Use \code{\link{posterior}} instead.
  
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
\keyword{classes}


