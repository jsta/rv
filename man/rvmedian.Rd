\name{rvmedian}
\title{Medians of the Components of a Random Vector}
\alias{rvmedian}
\usage{
  rvmedian(x, \dots)
}
\description{
  \code{rvmedian}
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{further arguments passed to \code{median}}
}
\details{
  \code{rvmedian} applies \code{median} to the simulations of
  each component of \code{x}.  
}
\value{
  A \emph{numeric} vector of the same dimension as \code{x}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\seealso{
  \code{\link{rvmin}}, \code{\link{rvmax}}, \code{\link{rvmean}}.
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

