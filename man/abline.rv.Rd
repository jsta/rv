\name{abline}
\title{Add (Random) Straight Lines to a Plot}
\alias{abline}
\usage{
abline (a = NULL, b = NULL, h = NULL, v = NULL, \dots)
}
\description{
  \code{abline}, with random arguments,
  plots a sample of lines corresponding to 
  of simulations of rv object \code{x}.
}
\arguments{
  \item{a}{...}
  \item{b}{...}
  \item{h}{...}
  \item{v}{..}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  This is a version of \code{abline} that accepts random variable objects
  for the arguments \code{a}, \code{b}, \code{h}, or \code{v}. 

  The number of lines is determined by 
  \code{rvpar("line.sample")}, default 20.

  The original help page is here: \code{\link[pkg:graphics]{abline}}.
 
  EXPERIMENTAL.
 
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
     demo("rvexample1")
  }
}

\keyword{aplot}

