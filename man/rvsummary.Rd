\name{rvsummary}
\title{Random Vector Summaries}
\alias{is.rvsummary}
\alias{as.rvsummary}
\alias{as.rvsummary.rv}
\alias{as.rvsummary.rvsummary}
\alias{print.rvsummary}
\alias{print.rvsummary_rvfactor}
\alias{dim.rvsummary}
\alias{dim<-.rvsummary}
\alias{dimnames.rvsummary}
\alias{dimnames<-.rvsummary}
\alias{names.rvsummary}
\alias{names<-.rvsummary}
\alias{length.rvsummary}
\alias{as.double.rvsummary}
\usage{
  is.rvsummary(x)
  as.rvsummary(x, \dots)
  \method{as.rvsummary}{rv}(x, quantiles = (0:200/200), \dots)
  \method{as.rvsummary}{rvsummary}(x, \dots)
  \method{print}{rvsummary_rvfactor}(x, all.levels=FALSE, \dots)
  \method{dim}{rvsummary}(x)
  \method{dim}{rvsummary}(x) <- value
  \method{dimnames}{rvsummary}(x)
  \method{dimnames}{rvsummary}(x) <- value
  \method{names}{rvsummary}(x)
  \method{names}{rvsummary}(x) <- value
  \method{length}{rvsummary}(x)
  \method{as.double}{rvsummary}(x, \dots)
}
\description{
  \code{rvsummary} is a class of objects
  that hold the summary information on each scalar component of a random variable
  (quantiles, mean, sd, number of simulations etc.)
}
\arguments{
  \item{x}{object to be coerced or tested}
  \item{quantiles}{quantiles to calculate and store in the object}
  \item{all.levels}{logical; whether to print all levels or not (see below for details)}
  \item{value}{either \code{NULL} or a numeric vector whose product matches the length of the \code{rvsummary} object}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The \code{rvsummary} class provides a means to store a concise representation
  of the marginal posterior distributions of the vector components.
  By default, the 201 quantiles
  \preformatted{
     0, 0.005, 0.01, 0.015, ..., 0.990, 0.995, 1 
  }
  are saved for each vector component in an \code{rvsummary} object.

  \code{is.rvsummary} tests whether the object is an \code{rvsummary} object;
  \code{as.rvsummary} coerces a random vector object to a \code{rvsummary} object.

  There following (generic) functions work with \code{rvsummary} objects:
  \code{rvmean},
  \code{rvsd},
  \code{rvvar},
  \code{rvquantile},
  \code{rnsims},
  \code{sims},
  and consequently any `rv-only' function that depends only on these functions
  will work; e.g. \code{is.constant}, which depends only on \code{rvnsims}.

  The method \code{is.double} is provided for compatibility reasons;
  this is needed in a function called by \code{plot.rvsummary}

  The arithmetic operators and mathematical functions will not work with
  \code{rvsummary} objects. 

  The \code{sims} method returns the quantiles. 
}
\value{
  An object of class \code{rvsummary} \emph{and} of subclass
  \code{rvsummary_numeric},
  \code{rvsummary_integer},
  \code{rvsummary_logical},
  or
  \code{rvsummary_rvfactor}.
}
\note{

}
\seealso{
  \code{\link{rvfactor}}
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

\examples{
  x <- rvnorm(mean=1:12)
  sx <- as.rvsummary(x)
  print(sx)          # prints the summary of the rvsummary object
  length(sx)         # 12
  dim(sx)            # NULL
  dim(sx) <- c(3,4)  #   
  dimnames(sx) <- list(1:3, 1:4)
  names(sx) <- 1:12  # 
  print(sx)          # prints the names and dimnames as well  
}
\keyword{classes}

