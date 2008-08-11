\name{rvquantile}
\title{Componentwise Quantiles of Random Variables}
\alias{rvquantile}
\alias{rvquantile.rv}
\alias{rvquantile.rvsummary}
\alias{rvmedian}
\usage{
  rvquantile(x, \dots)
  \method{rvquantile}{rv}(x, probs=c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975), ignoreInf=FALSE, \dots)
  \method{rvquantile}{rvsummary}(x, probs=c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975), \dots)
  rvmedian(x)
}
\description{
  Computes componentwise quantiles of random vectors or arrays.
}
\arguments{
  \item{x}{an object}
  \item{probs}{numeric vector of probabilities with values in \emph{[0,1]}}
  \item{ignoreInf}{ignore infinite values}
  \item{\dots}{further arguments passed to \code{quantile}}
}
\details{
  \code{rvquantile} applies the \code{quantile} function to each
  column of \code{sims(x)}.

  \code{rvmedian} applies \code{median} to the each column of \code{sims(x)}.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\value{
  A \emph{numeric} vector of quantiles.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  x <- rvnorm(3)
  rvquantile(x)
  rvquantile(x, probs=c(0, 0.01, 0.99, 1))
  rvmedian(x)
}

\keyword{classes}

