\name{rvlapply}
\title{Apply a Function to Components of a Random Vector}
\alias{lapply.rv}
\usage{
  lapply.rv(x, FUN, \dots)
}
\description{
  \code{lapply.rv} applies a given function to 
  the components of a random vector; it is an `rv'-compatible
  version of \code{lapply}.
}
\arguments{
  \item{x}{an object}
  \item{FUN}{a function}
  \item{\dots}{further arguments passed to the function \code{FUN}}
}
\details{

}
\value{
  A list of the same length as \code{x}, with values
  \code{FUN(x[i], ...)} in component \code{i}.
  The names of the list are taken from the names of the components of \code{x}.
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

