\name{as.bugs.rv}
\title{Coerce an rv object to a bugs object}
\alias{as.bugs.rv}
\usage{
  as.bugs.rv (x, DIC = FALSE)
}
\description{
  Coerces an rv object to a bugs object.
}
\arguments{
  \item{x}{an rv}
  \item{DIC}{logical, if TRUE, compute deviance information criterion (DIC)}
}
\value{
  a bugs object.
}
\note{
  Works currently ONLY with random variable vectors that have
  been created from 3-dimensional arrays containing the MCMC chains.
  This may change in the future.
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
%\examples{
% #
%}

\keyword{classes}

