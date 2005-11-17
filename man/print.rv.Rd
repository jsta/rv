\name{print.rv}
\title{Print Distribution Summaries of Random Variables}
\alias{print.rv}
\alias{summary.rv}
\alias{print.rvsim}
\usage{
print.rv(x, \dots)
}
\description{
  \code{print.rv}
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  \code{print.rv} 
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}

\examples{
  print(rvnorm(mean=rvnorm(1)))
}
\keyword{classes}

%%%%%
