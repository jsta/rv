\name{Extremes-rv}
\title{Maxima and Minima of Random Variables}
\alias{rvmin}
\alias{rvmax}
\alias{rvrange}
\alias{min.rv}
\alias{max.rv}
\alias{pmin.rv}
\alias{pmax.rv}
\usage{
  rvmin(x)
  rvmax(x)
  rvrange(x)
  \method{min}{rv}(\dots, na.rm=FALSE)
  \method{max}{rv}(\dots, na.rm=FALSE)
  \method{pmin}{rv}(\dots, na.rm=FALSE)
  \method{pmax}{rv}(\dots, na.rm=FALSE)
}
\description{
  Returns the maxima and minima of the components of a random vector.
}
\arguments{
  \item{x}{an \code{rv} or \code{rvsummary} object}
  \item{na.rm}{remove missing values?}
  \item{\dots}{one or more \code{rv} objects or numeric objects}
}
\details{
  \code{rvmin} applies the function \code{min} to each component of the argument \code{x}.
  Missing values are removed. 

  \code{rvmax} applies the function \code{max} to each component of the argument \code{x}.
  Missing values are removed. 

  \code{rvrange} applies the function \code{range} to each component of
  the argument \code{x}.  Missing values are removed. 

  \code{min.rv} returns the minimum of the random \emph{vector},
  returning thus one random variable. Similarly \code{max.rv} returns
  the maximum of a vector.

  \code{pmin.rv} and \code{pmax.rv} returns the componentwise minima or
  maxima of several random vectors or constants, yielding thus a random
  vector of the same length. 

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

