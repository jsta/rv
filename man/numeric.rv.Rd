\name{numeric.rv}
\title{Numeric Random Vectors}
\alias{numeric.rv}
\alias{as.numeric.rv}
\alias{is.numeric.rv}
\alias{is.numeric.rvfactor}
\alias{as.numeric.rvfactor}
\usage{
  \method{is.numeric}{rv}(x)
  \method{as.numeric}{rv}(x, \dots)
  \method{is.numeric}{rvfactor}(x)
  \method{as.numeric}{rvfactor}(x, \dots)
}
\description{
  Creates or coerces rv objects of type "numeric". 
}
\arguments{
  \item{x}{an rv object to be coerced or tested.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  \code{is.numeric(x)} returns \code{TRUE} if and only if \emph{each} component 
  of \code{x} is numeric-valued.

   \code{as.numeric.rv} coerces an rv object into numeric-valued one.
   In effect, the function  \code{as.numeric} is applied to all simulations.

  Random factors are not numeric (just as non-random factors aren't).
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\note{

}
\seealso{
  \code{\link{numeric}}.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  x <- as.logical(rvbern(1,0.5)) # Bernoulli rv
  is.numeric(x)           # FALSE
  x <- as.numeric(x)      # coerce to numeric; all TRUEs become ones, FALSEs zeros
  is.numeric(x)           # TRUE
}

\keyword{classes}

