\name{points.rv}
\title{Add Points and Intervals to a Plot}
\alias{points.rv}
\usage{
points.rv(x, y = NULL, what=NULL, rvcol=c("green", "red"), \dots)
}
\description{
  Draw a sequence of points or uncertainty intervals at specified
  (fixed) x-coordinates.
}
\arguments{
  \item{x}{x-coordinates}
  \item{y}{y-coordinates}
  \item{what}{what to plot: intervals, means, medians...}
  \item{rvcol}{colors for the intervals}
  \item{\dots}{further arguments passed to \code{points}}
}
\details{
  (Under development.)
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{classes}


