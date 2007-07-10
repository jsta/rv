\name{print.rv}
\title{Print Distribution Summary of a Random Variable}
\alias{print.rv}
\alias{print.rvfactor}
\alias{print.rvmixed}
\alias{print.rvsim}
\usage{
  \method{print}{rv}(x, digits=NULL, \dots)
  \method{print}{rvfactor}(x, digits=NULL, \dots)
  \method{print}{rvmixed}(x, digits=NULL, \dots)
}
\description{
  Prints a summary of the random variable object.
}
\arguments{
  \item{x}{an rv object}
  \item{digits}{minimal number of significant digits, default \code{rvpar('print.digits')}}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  Invokes first the summary method of the object, then prints the result.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\seealso{
  \code{\link{summary.rv}}
}
\examples{
  print(rvnorm(mean=rvnorm(1)))
}
\keyword{classes}

%%%%%
