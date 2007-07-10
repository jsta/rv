\name{rvinci}
\title{Is a value inside a confidence interval?}
\alias{rvinci}
\usage{
  rvinci(obj, x, interval=0.95, one.sided=FALSE, left=FALSE)
}
\description{
  \code{rvinci} returns \code{TRUE} iff a given value is inside
  a given interval.
}
\arguments{
  \item{obj}{a random variable.}
\item{x}{a numeric object.}
\item{interval}{a number between 0 and 1.}
\item{one.sided}{logical. Is the interval one-sided?}
\item{left}{logical. Is the one-sided interval on the left?}
}
\details{
  Returns \code{TRUE} if \code{x} (a number) is in the \code{100*interval}\% interval 
  of the random variable \code{obj}.
}
\value{
  A logical vector.
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
  #
}

\keyword{classes}

