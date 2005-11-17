\name{hist.rv}
\title{Histogram of a random vector}
\alias{hist.rv}
\usage{
 \method{hist}{rv}(x, grid = c(4, 5), xlim = x.range, main = paste(xname, "simulation"), \dots)
}
\description{
 \code{hist.rv} shows a grid of histograms generated
 from random draws of the random vector argument.
}
\arguments{
  \item{x}{an object}
\item{grid}{a vector of two numbers, indicating the size of the grid to plot the histograms}
\item{xlim}{x limits}
\item{main}{main title}
\item{\dots}{Other arguments passed on to \link{hist}}
}
\details{
  
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  \dontshow{require(rv)}
  x <- rvnorm(30)
  hist(x)
}

\keyword{classes}

