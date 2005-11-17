\name{rowapply}
\title{Apply a Function to Rows of a Random Matrix}
\alias{rowapply}
\usage{
  rowapply(X, FUN, \dots)
}
\description{
  \code{rowapply} takes a random matrix and applies a function 
  to its rows, returning another rv object.
}
\arguments{
  \item{X}{random matrix (an rv object with non-null \code{dim}ension)}
  \item{FUN}{a function to apply}
  \item{\dots}{further arguments passed to \code{FUN}}
}
\value{
  A random variable (an rv object).
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


