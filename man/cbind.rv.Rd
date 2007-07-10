\name{cbind.rv}
\title{Combine random vectors by columns or rows}
\alias{cbind.rv}
\alias{rbind.rv}
\usage{
 \method{cbind}{rv}(\dots, deparse.level = 1)
 \method{rbind}{rv}(\dots, deparse.level = 1)
}
\description{
  Combines random vectors by columns (\code{cbind.rv}) or rows (\code{rbind.rv}).
}
\arguments{
  \item{\dots}{vectors or matrices, can be rv objects}
  \item{deparse.level}{(passed on to cbind)}
}
\details{
  See \link{cbind} and \link{rbind} for details.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvnorm(10)
  y <- rvnorm(10)
  cbind(x,y)
  rbind(x,y)
}

\keyword{classes}

