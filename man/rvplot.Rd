\name{rvplot}
\title{Scatterplot of Simulations of a Random Vector}
\alias{rvplot}
\usage{
rvplot(x,y,\dots)
}
\description{
  Produces a 2-dimensional scatterplot of the simulations of a random vector. 
}
\arguments{
  \item{x}{an rv object; only two first components are used}
  \item{y}{(optional) rv object}
  \item{\dots}{further arguments passed to \code{plot}}
}
\details{
  Only the two first scalars are used if only \code{x} is specified.
  If both \code{x} and \code{y} are specified, they are assumed to be scalars.
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

