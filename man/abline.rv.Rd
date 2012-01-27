\name{abline.rv}
\title{Add (Random) Straight Lines to a Plot}
\alias{abline.rv}
\usage{
  abline.rv (a = NULL, b = NULL, h = NULL, v = NULL, \dots)
}
\description{
  \code{abline.rv}, with random arguments (i.e. arguments of which at
  least one is an \code{rv} object),
  plots a sample of lines corresponding to 
  of simulations of rv object \code{x}.
  If the arguments are all numeric (none is an \code{rv} object),
  the function call is passed on to \code{abline}.
}
\arguments{
  \item{a}{intercept}
  \item{b}{slope}
  \item{h}{y-value(s) horizontal line(s)}
  \item{v}{x-value(s) horizontal line(s)}
  \item{\dots}{further arguments passed to \code{\link{abline}}}
}
\details{
  This is a version of \code{abline} that accepts random variable objects
  for the arguments \code{a}, \code{b}, \code{h}, or \code{v}. 

  The number of lines is determined by 
  \code{rvpar("line.sample")}, default 20.

  See the original help page in package `graphics.'
 
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

