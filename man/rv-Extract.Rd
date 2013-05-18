\name{Extract.rv}
\title{Extract or Replace Parts of a Random Vector}
\alias{[.rv}
\alias{[.rvfactor}
\alias{[.rvsummary}
\alias{[<-.rv}
\alias{[<-.rv}
\alias{[<-.rvsummary}
\alias{impute<-}
\usage{
\method{[}{rv}(x, \dots , drop = TRUE)
\method{[}{rvfactor}(x, \dots , drop = FALSE)
\method{[}{rvsummary}(x, \dots , drop = TRUE)
\method{[}{rv}(x, \dots) <- value
\method{[}{rvsummary}(x, \dots) <- value
impute(x, \dots) <- value
}
\arguments{
  \item{x}{object from which to extract element(s) or in which to
    replace element(s).}
  \item{\dots}{indices specifying elements to extract or replace.}
  \item{value}{typically an array-like R object of a similar class as \code{x}.}
  \item{drop}{For matrices and arrays.  If \code{TRUE} the result is 
    coerced to  the lowest possible dimension (see the examples).
    This only works for extracting elements, not for the replacement.}
}
\description{
  Bracket slice and assignment methods adapted for random vectors and
arrays.
  The assignment function \code{impute<-} is compatible with both non-rv
and rv objects (rv, rvsummary, and rvfactor objects). To write universal
code that works both atomic and rv objects, use \code{impute(x, ...) <-
value} instead of \code{x[...] <- value}.
}
\details{
  
  NOTE. \code{x} will NOT be automatically coerced
  into an rv object.
  
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
%\seealso{
%  See \link{rvcompatibility} for examples.
%}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}

\examples{
  x <- rvnorm(1)
  y <- (1:5)
  \dontrun{
     y[2] <- x ## Will not work
  }
  impute(y, 2) <- x
}
\keyword{classes}

