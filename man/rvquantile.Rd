\name{rvquantile}
\title{Componentwise Quantiles of Random Variables}
\alias{rvquantile}
\usage{
  rvquantile(x, \dots)
}
\description{
  Computes componentwise quantiles of random vectors or arrays.
}
\arguments{
  \item{x}{an object}
  \item{\dots}{further arguments passed to \code{quantile}}
}
\details{
  \code{rvquantile} applies the \code{quantile} function to each
  column of \code{sims(x)}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\value{
  A \emph{numeric} vector of quantiles.
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

