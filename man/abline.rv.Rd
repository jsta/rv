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
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  \dontrun{
     demo("rvexample1")
  }
}

\keyword{aplot}

