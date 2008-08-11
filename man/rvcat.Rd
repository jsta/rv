\name{rvcat}
\title{Generate Categorical Random Variables}
\alias{rvcat}
\usage{
  rvcat(n=1, prob, levels=NULL)
}
\description{
  Generates a random factor (i.e. a categorical random variable),
  given the probabilities of each category and their corresponding labels.
}
\arguments{
  \item{n}{integer, number of random variables to generate}
  \item{prob}{vector of probabilities of successes of each trial (may be constant or an rv object)}
  \item{levels}{(character) labels for the categories}
}
\details{
  The length of \code{prob} determines the number of bins.

  The vector \code{prob} will be normalized to have sum 1.

}
\value{
  A \emph{random factor} of length \code{length(prob)}. 
}
\seealso{
  \code{\link{rvfactor}}
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
  rvcat(1, prob=c(0.5, 0.3, 0.2)) # default levels: 1, 2, 3
  rvcat(1, prob=c(5, 3, 2)) # same as above
  p <- rvdirichlet(1, alpha=c(0.7, 0.3)) # prior probabilities
  rvcat(1, prob=p, levels=c("Group 1", "Group 2"))
}

\keyword{classes}

