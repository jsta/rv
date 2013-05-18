\name{hist.rv}
\title{Histogram of a random vector}
\alias{hist.rv}
\usage{
 \method{hist}{rv}(x, grid = c(4, 5), xlim = x.range,
main = paste(xname, "simulation"), freq=FALSE, \dots)
}
\description{
 \code{hist.rv} shows a grid of histograms generated
 from random draws of the random vector argument.
}
\arguments{
  \item{x}{an object}
  \item{grid}{a vector of two numbers, indicating the size of the grid
    to plot the histograms}
\item{xlim}{x limits}
\item{main}{main title}
\item{freq}{logical; if \code{FALSE}, plots as probability density,
  as it should.}
\item{\dots}{Other arguments passed on to \link{hist}}
}
%\details{
%  
%}
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
  \dontshow{require(rv)}
  x <- rvnorm(30)
  hist(x)
}

\keyword{classes}

