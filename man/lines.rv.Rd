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
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{ 
  x <- 1:10
  y <- rvnorm(mean=x)
  par(mfrow=c(2,2))
  plot(x, y, type="b", main="Intervals and random lines", rvcol="blue", col="gray")
  plot(x, y, type="l", main="Only random lines", col="gray")
  plot(as.rv(x), E(y), type="b", main="Means, connected by a constant line",
     rvcol="blue", col="gray")
  plot(y, x, type="l", main="Only random lines", col="gray")
}

\keyword{aplot}

%%% 
