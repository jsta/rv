\name{print.rv}
\title{Print Distribution Summary of a Random Variable}
\alias{print.rv}
\usage{
  \method{print}{rv}(x, digits=rvpar("print.digits"), \dots)
}
\description{
  Prints a summary of the random variable object.
}
\arguments{
  \item{x}{an rv object}
  \item{digits}{minimal number of significant digits}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  Invokes first the summary method of the object, then prints the result.
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
  \code{\link{summary.rv}}, \code{\link{rvfactor}}
}
\examples{
  print(rvnorm(mean=rvnorm(1)))
}
\keyword{classes}

%%%%%
