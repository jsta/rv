% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuzzy.R
\name{is.fuzzy}
\alias{is.fuzzy}
\alias{fuzzy}
\alias{is.fuzzy.rv}
\title{Fuzziness}
\usage{
is.fuzzy(x)
}
\arguments{
\item{x}{an object, random or constant}
}
\description{
Tests whether an object is "fuzzy", i.e.  a logical random scalar that has
probability strictly between zero and one (not strictly true nor strictly
false).
}
\examples{

  x <- as.logical(rvbern(1,0.4)) # a logical random variable
  is.fuzzy(x) # TRUE, since x is logical and not constant
  is.fuzzy(x<2) # FALSE, since x is less than 2 with probability one
  is.fuzzy(rvnorm(1)) # FALSE, since it's not a probability
  is.fuzzy(TRUE) # FALSE, since TRUE is strictly TRUE
  is.fuzzy(1) # FALSE, since 1 is not a logical variable

}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects. Statistics and
Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\author{
Jouni Kerman \email{jouni@kerman.com}
}
\keyword{classes}
