\name{rvmean}
\title{Expectation of a Random Variable}
\alias{rvmean}
\alias{E}
\alias{Pr}
\usage{
  rvmean(x)
  E(x)
  Pr(X)
}
\description{
  \code{rvmean}
}
\arguments{
  \item{x}{an rv object}
  \item{X}{a logical rv object}
}
\details{
\code{rvmean} computes the means of the simulations 
of all individual components of a random vector (rv) object.

\code{E} is an alias for \code{rvmean}, standing for ``Expectation.''

\code{Pr} is another alias for \code{rvmean}, standing for ``Probability of'';
suggested to be used when the argument is a logical statement involving random variables
(that is, a description of an event such as \code{x>0} or \code{x>y}).
Then \code{Pr(x>0)} gives the probability of the event ``x>0''.
The statement \code{x>0} returns a Bernoulli
(indicator) random variable object (having 1/0 or TRUE/FALSE values)
and the expectation of such variable is just the probability of the 
event where the indicator is one.
}
\value{
  A \emph{numerical} vector with the same dimension as \code{x}.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\seealso{
  \code{\link{mean.rv}}: distribution of the arithmetic mean of a vector;
  \code{\link{rvmin}}, \code{\link{rvmax}}, \code{\link{rvmedian}}, 
  \code{link{rvvar}}, \code{\link{rvsd}}.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  x <- rvnorm(mean=(1:10)/5, sd=1)
  rvmean(x)  # means of the 10 components
  E(x)       # same as rvmean(x)
  Pr(x>1)    # probabilities that each component is >1.
}
\keyword{classes}

