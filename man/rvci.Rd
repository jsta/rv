\name{rvci}
\title{Credible (Uncertainty) Intervals for Random Scalars}
\alias{rvci}
\usage{
  rvci(obj, interval=0.95, one.sided=FALSE, left=TRUE)
}
\description{
  Computes credible (uncertainty) intervals for a given vector.
}
\arguments{
  \item{obj}{random scalar or vector}
  \item{interval}{a number between 0 and 1}
  \item{one.sided}{logical, FALSE if two-sided interval is desired}
  \item{left}{logical, indicating if the left one-sided interval is desired}
}
\value{
  For two-sided intervals,
  an array of numbers of dimension \code{c(2,length(x))},
  for one-sided intervals,
  a vector of the same length as \code{x}.
}
\references{ 
Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations 
Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  rvci(rvnorm(1), interval=0.683) # Should be about c(-1,1).
}

\keyword{classes}

