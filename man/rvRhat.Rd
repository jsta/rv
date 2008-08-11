\name{rvRhat}
\title{R-hat Convergence Diagnostic}
\alias{rvRhat}
\usage{
rvRhat(x)
}
\description{
  Retrieves the R-hat convergence diagnostic for each component of the argument
}
\arguments{
  \item{x}{an object}
}
\details{
  The R-hat values are assumed to be saved as attributes.
  If they are not available, \code{NA} will be returned.

  R-hat is computed by programs such as Umacs and R2WinBUGS.
}
\value{
  Vector of numbers, \code{NA} if R-hat is not available.
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

\keyword{classes}

