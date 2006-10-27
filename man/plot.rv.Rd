\name{plot.rv}
\title{Plotting Random Vectors}
\alias{plot.rv}
\usage{
  \method{plot}{rv}(x, y, what=c("95\%","50\%","mean","median"), ylim=range(sims(y)), xlim=range(sims(x)), rvcol=c("grey20", "grey40"), \dots)
}
\description{
  Draw a "random scatter plot". 
  Currently only the y coordinate can be random,
  so intervals are plotted for each fixed x-coordinate..
}
\arguments{
  \item{x}{(optional) fixed x-coordinate}
  \item{y}{random vector}
  \item{what}{what to plot: an interval or a dot}
  \item{xlim}{the x limits of the plot}
  \item{ylim}{the y limits of the plot}
  \item{rvcol}{colors for the intervals (passed on to \code{link{points.rv}})}
  \item{\dots}{other arguments passed on to \code{plot}}
}
\details{
  (This function is still under development.)
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- 1:30
  y <- rvnorm(mean=x, sd=1)
  \dontrun{plot.rv(x,y)}  # x is not an rv, so be explicit and specify plot.rv 
  \dontrun{plot(y)}       # y is an rv, so R knows to call plot.rv.
}

\keyword{classes}

