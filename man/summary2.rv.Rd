\name{summary.rv}
\title{Summary of a Random Variable}
\alias{summary.rv}
\alias{summary.rvfactor}
\alias{summary.rvmixed}
\usage{
  \method{summary}{rv}(object, \dots)
  \method{summary}{rvfactor}(object, all.levels=FALSE, \dots)
  \method{summary}{rvmixed}(object, all.levels=FALSE, \dots)
}
\description{
  Gives a numerical summary of the random variable. 
}
\arguments{
  \item{object}{object to summarize}
  \item{all.levels}{show summary for all levels even if there are too many to display in one line}
  \item{\dots}{rv object}
}
\details{
      
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvarray(rvnorm(6), c(2,3))
  summary(x)
  summary(rvfactor(trunc(x)))
}

\keyword{internal}

