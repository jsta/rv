% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvexp.R
\name{rvexp}
\alias{rvexp}
\title{Generate Random Vectors from an Exponential Sampling Model}
\usage{
rvexp(n = 1, rate = 1)
}
\arguments{
\item{n}{integer: number of variables to generate}

\item{rate}{prior distribution for the rate parameter (constant or random)}
}
\description{
\code{rvexp}
}
\details{
\code{rvexp}
}
\examples{

  y <- rvexp(1, rate=rvexp(1)) # What marginal distribution does y have now?

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
