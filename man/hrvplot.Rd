\name{hrvplot}
\title{Horizontal interval plot of components of a random vector}
\alias{hrvplot}
\alias{hrvplot.matrix}
\alias{hrvplot.rv}
\usage{
hrvplot(x, labels=NULL, pos=NULL, output.file=NULL)
\method{hrvplot}{matrix}(x, labels, pos=NULL, output.file=NULL)
\method{hrvplot}{rv}(x, labels=NULL, \dots)
}
\description{
  \code{hrvplot} plots uncertainty (that is, ``credible") intervals 
  for the scalar components of the given random vector.
}
\arguments{
  \item{x}{an object}
\item{labels}{alternative labeling}
\item{pos}{}
\item{output.file}{alternative output file name}
\item{\dots}{}
}
\details{
  \code{hrvplot} 
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

