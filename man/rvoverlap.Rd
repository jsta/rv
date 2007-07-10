\name{rvoverlap}
\title{Test Whether Intervals Overlap with a Given Interval}
\alias{rvoverlap}
\usage{
  rvoverlap(obj, iv, interval=0.95, one.sided=FALSE, left=TRUE)
}
\description{
  Tests whether the uncertainty intervals of
  the given random variables overlap with a given interval.

  NOTE. Still experimental.
}
\arguments{
  \item{obj}{random variable (vector or array)}
  \item{iv}{an interval (range)}
  \item{interval}{size of the (middle) interval }
  \item{one.sided}{logical, FALSE if two-sided interval is desired}
  \item{left}{logical, indicating if the left one-sided interval is desired}
}
\value{
  \item{0}{Middle \code{interval} is fully enclosed by \code{iv}}
  \item{-1}{Middle \code{interval} is partially covered on the left of \code{iv}}
  \item{+1}{Middle \code{interval} is partially covered on the right of \code{iv}}
  \item{-Inf}{Middle \code{interval} is completely to the left of \code{iv}}
  \item{+Inf}{Middle \code{interval} is completely to the right of \code{iv}}
}
\references{ 
Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations 
Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  
  rvoverlap(rvnorm(20), iv=c(-1,1), interval=0.683) # 
}

\keyword{classes}

