\name{detachrv}
\title{Detach the rv package}
\alias{detachrv}
\usage{
  detachrv()
}
\description{
  \code{detachrv} detaches the rv package and restores the original
  functions in \code{base}, \code{graphics} and \code{stats} packages.
}
%\arguments{
%  
%}
\details{
  Currently equivalent to \code{detach("package:rv")}.
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
  \dontrun{
   library(rv)
   detachrv()
  }
}

\keyword{classes}


