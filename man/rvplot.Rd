\name{rvplot}
\title{Horizontal interval plot of components of a random vector}
\alias{rvplot}
\alias{rvplot.matrix}
\alias{rvplot.rv}
\usage{
rvplot(x, labels=NULL, pos=NULL, output.file=NULL, \dots)
\method{rvplot}{matrix}(x, labels, pos=NULL, output.file=NULL, xlim=NULL, vline=0, \dots)
\method{rvplot}{rv}(x, labels=NULL, \dots)
}
\description{
  \code{rvplot} plots uncertainty (``credible") intervals 
  for the scalar components of the given random vector.
}
\arguments{
  \item{x}{an object}
\item{labels}{alternative labeling}
\item{pos}{}
\item{output.file}{alternative output file name}
\item{xlim}{x-limits}
\item{vline}{x-coordinate of the vertical reference line}
\item{\dots}{}
}
\note{
  Under development!
}
%\details{
%  \code{rvplot} 
%}
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

