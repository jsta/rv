\name{lines.rv}
\title{Add Connected (Random) Line Segments to a Plot}
\alias{lines.rv}
\usage{
  \method{lines}{rv}(x, y, type="l", \dots)
}
\description{
  Adds a sample of line segments randomly drawn from the joint distribution of 
  \code{(x,y)}. 
}
\arguments{
  \item{x, y}{coordinate vectors of points to join}
  \item{type}{character indicating the type of plotting, currently 'l' and 'p' are the only possibilities}
  \item{\dots}{further arguments passed to \code{points}}
}
\details{
  The size of the sample (number of segments drawn) 
  is determined by
  \code{rvpar(line.sample)}.

  \code{lines.rv} is implemented as part of \code{points.rv}. 

  See \code{\link{points.rv}} for details of the parameters.
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
  x <- as.rv(1:10)
  y <- rvnorm(mean=x)
  par(mfrow=c(2,2))
  plot(x, y, type="b", main="Intervals and random lines", rvcol="blue", col="gray")
  plot(x, y, type="l", main="Only random lines", col="gray")
  plot(x, E(y), type="b", main="Means, connected by a constant line", col="gray")
  plot(x, rvmedian(y), type="b", pch=19, main="Median & middle 95 pc CI band", col="darkgray")
  lines(rvquantile(y, 0.025), col="gray")
  lines(rvquantile(y, 1-0.025), col="gray")
}

\keyword{aplot}

%%% 
