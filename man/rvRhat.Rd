\name{rvRhat}
\title{R-hat Convergence Diagnostic}
\alias{rvRhat}
\usage{
rvRhat(x)
}
\description{
  Retrieves the R-hat convergence diagnostic for each component of the argument
}
\arguments{
  \item{x}{an object}
}
\details{
  The R-hat values are assumed to be saved as attributes.
  If they are not available, \code{NA} will be returned.
}
\value{
  Vector of numbers, \code{NA} if R-hat is not available.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}

\keyword{classes}

