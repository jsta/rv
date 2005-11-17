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
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvnorm(mean=1:3)
  y <- rvnorm(mean=2:4)
  rvcov(x,y)
  rvcov(x,x)
}

\keyword{classes}

