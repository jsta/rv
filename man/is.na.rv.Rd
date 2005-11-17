\name{is.na.rv}
\title{Distribution of the Indicator Function of Missing Data}
\alias{is.na.rv}
\usage{
  \method{is.na}{rv}(x)
}
\description{
  \code{is.na.rv}
}
\arguments{
  \item{x}{an rv object}
}
\details{
  Internally, \code{is.na.rv} applies the function \code{is.na}
  to each simulation of each component of the argument vector.
}
\value{
  A ``Bernoulli" random vector of the same length and dimension as those of \code{x}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
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

%%% maybe make up an example where we have random-length Poisson variable?
