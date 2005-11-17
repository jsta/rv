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
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  y <- rnorm(30) # Some data.
  x <- rvboot(y) # A random vector of length 30 (each component has the same distribution)
  print(mean(x)) # Bootstrap estimate of the mean.
  print(sd(x))   # Bootstrap estimate of the sd.
  rvinci(mean(x), 0) # Hypothesis test: mean of x is zero (at 5% c.l.). FALSE => reject.
}

\keyword{classes}

