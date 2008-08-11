\name{rvboot}
\title{Generate a Random Vector from an Empirical Distribution}
\alias{rvboot}
\usage{
  rvboot(data)
}
\description{
  \code{rvboot} generates a random vector of the same length as data
  from the empirical distribution of the data.
}
\arguments{
  \item{data}{A vector of constants}
}
\details{
  \code{rvboot} 
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
  y <- rnorm(30) # Some data: 30 draws from standard normal.
  x <- rvboot(y) # A random vector of length 30 (each component has the same distribution)
  print(mean(x)) # Bootstrap estimate of the mean.
  print(sd(x))   # Bootstrap estimate of the sd.
%  rvinci(mean(x), 0) # Hypothesis test: mean of x is zero (at 5% level) FALSE => reject.
}

\keyword{classes}

