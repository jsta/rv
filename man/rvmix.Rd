\name{rvmix}
\title{Generate Random Variables from a Mixture Distribution}
\alias{rvmix}
\usage{
  rvmix(n=1, prob=NULL, index=NULL, components=list())
}
\description{
  Generates a random variable 
  from a mixture of given random variables.
}
\arguments{
  \item{n}{number of independent draws to output}
  \item{prob}{probabilities of the components (or weights, will be normalized)}
  \item{index}{(optional) index random variable, giving the distribution of indices}
  \item{components}{a list of the component random variables}
}
\details{
  If an index is given, \code{prob} is ignored.
}
\value{
  A random variable.
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
  comp <- list(rvnorm(mean=0), rvnorm(mean=10)) # Normal components of a mixture distribution
  rvmix(prob=c(95,5), components=comp) # 5 per cent comes from N(10,1)
}

\keyword{classes}

