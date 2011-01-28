\name{Extract.rv}
\title{Extract or Replace Parts of a Random Vector}
\alias{[.rv}
\alias{[.rvsummary}
\alias{[<-.rv}
\alias{[<-.rv}
\alias{[<-.rvsummary}
\synopsis{
%x[i, j, \dots, drop = TRUE]
%x[i, j, \dots , drop = TRUE] <- value
%x[[i, j, \dots, drop = TRUE]]
%x[[i, j, \dots , drop = TRUE]] <- value
\method{[}{rv}(x, \dots , drop = TRUE)
\method{[}{rvsummary}(x, \dots , drop = TRUE)
\method{[<-}{rv}(x, \dots , drop = TRUE)
\method{[<-}{rvsummary}(x, \dots , drop = TRUE)
%\method{[[<-}{rv}(x, \dots , drop = TRUE)
}
\description{
  Bracket slice and assignment methods adapted for random vectors and arrays.
}
\details{
  
  \code{x} may also be a regular vector that is automatically coerced
  into an rv object if the compatibility flag has been 
  set to 1 by \code{rvcompatibility(1)} (this is the default).
  
  \code{value} may be an rv object or a regular numeric object.

  Extracting rv objects works the same way as extracting
  components of a numerical vector or array.
  The return value is always an object of class 'rv'.
  Type ?Extract for details.

  Note: the index arguments (\code{i}, \code{j}, etc.) 
  \emph{must} be constants, but this may change in the future.

  %Note: the index arguments (\code{i}, \code{j}, etc.) may be 
  %themselves random variables, however they will be coerced
  %into \emph{integers}, as one would expect.
}
\value{
 A random variable (an rv object).
}
\seealso{
  See \link{rvcompatibility} for examples.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}

\examples{
  #
}
\keyword{classes}

