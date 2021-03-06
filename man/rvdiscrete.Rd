% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvdiscrete.R
\name{rvdiscrete}
\alias{rvdiscrete}
\title{Generate Random Vectors from a Discrete Sampling Model}
\usage{
rvdiscrete(n = 1, x, prob = NULL)
}
\arguments{
\item{n}{integer: number of scalars to generate}

\item{x}{values of the distribution}

\item{prob}{probabilities (optional, default: all equal)}
}
\description{
Generates random variables from a discrete distribution (from a finite
population with replacement).
}
\details{
Computes a random vector of length \code{n}, consisting of identicallly
distributed discrete random scalars with the discrete distribution with
values \code{x} and corresponding probabilities \code{prob}.  If \code{prob}
is not given, all values are considered equally distributed.
}
\examples{
\dontrun{

  # 8 people draw a number each from 1..10 with replacement.
  # What is the probability that the highest number of the eight is "10"?
  u <- rvdiscrete(n=8, x=1:10) # 8 iid variables from the discrete uniform 1:10.
  Pr(max(u)==10)
  # What is the probability that the person with the 3rd smallest number
  # has at least "3"?
  s <- sort(u) # order distribution
  Pr(s[3]>=3)
  }

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
