\name{plot.rv}
\title{Plotting Scatterplots of Random Variable Objects}
\alias{plot.rv}
\usage{
  \method{plot}{rv}(x, y, \dots)
}
\description{
  Draw a "random scatter plot". 
}
\arguments{
  \item{x}{random or fixed vector}
  \item{y}{random or fixed vector}
  \item{\dots}{other arguments passed on to \code{plot}}
}
\details{
  If a component \code{x} is fixed and the corresponding component of \code{y} is random, 
  the resulting `point' is a vertical uncertainty ('credible') interval.

  If a component \code{y} is fixed and the corresponding component of \code{x} is random, 
  the resulting `point' is a horizontal uncertainty ('credible') interval.

  If a component of \code{x} and the corresponding component of \code{y} is random, 
  the resulting `point' is a scatterplot of simulations from the joint distribution of
  code{(x,y)}.

}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\seealso{
  \code{\link{mlplot}}
}
\examples{
  x <- 1:30
  y <- rvnorm(mean=x, sd=1)
  \dontrun{plot(x,y)}
  \dontrun{plot(y,x)}
  \dontrun{plot(y)}
}

\keyword{aplot}

