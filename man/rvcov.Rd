\name{rvcov}
\title{Covariance Between Components of Random Vectors}
\alias{rvcov}
\usage{
  rvcov(x, y=NULL, \dots)
}
\description{
  \code{rvcov}
}
\arguments{
  \item{x}{a random vector}
  \item{y}{(optional) a random vector}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  \code{rvcov} 
}
\value{
  A covariance matrix.
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
  x <- rvnorm(mean=1:3)
  y <- rvnorm(mean=2:4)
  rvcov(x,y)
  rvcov(x,x)
}

\keyword{classes}

