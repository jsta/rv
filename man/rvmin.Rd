\name{Extremes-rv}
\title{Maxima and Minima of Random Variables}
\alias{rvmin}
\alias{rvmax}
\alias{rvrange}
\usage{
  rvmin(x)
  rvmax(x)
  rvrange(x)
}
\description{
  Returns the maxima and minima of the components of a random vector.
}
\arguments{
  \item{x}{an \code{rv} or \code{rvsummary} object}
}
\details{
  \code{rvmin} applies the function \code{min} to each component of the argument \code{x}.
  Missing values are removed. 

  \code{rvmax} applies the function \code{max} to each component of the argument \code{x}.
  Missing values are removed. 

  \code{rvrange} applies the function \code{range} to each component of the argument \code{x}.  Missing values are removed. 
}
\value{
  A \emph{numeric} vector of the same dimension as \code{x}.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\seealso{
  \code{\link{rvmedian}}, \code{\link{rvmean}}.
}
\examples{
  x <- rvpois(10, lambda=3)
  rvmin(x)
  rvmax(x)
  rvrange(x)
}

\keyword{classes}

