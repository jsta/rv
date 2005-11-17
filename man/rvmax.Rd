\name{rvmax}
\title{Maximum Simulated Value of Random Variable}
\alias{rvmax}
\usage{
  rvmax(x, \dots)
}
\description{
  \code{rvmax}
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{further arguments passed to \code{max}}
}
\details{
  \code{rvmax} 
}
\value{
  A \emph{numeric} vector of the same dimension as \code{x}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\seealso{
  \code{\link{rvmin}}, \code{\link{rvmedian}}, \code{\link{rvmean}}.
}
\examples{
  rvmax(rvnorm(1))
}

\keyword{classes}

