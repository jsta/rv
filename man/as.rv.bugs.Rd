\name{as.rv.bugs}
\title{Coerce an R2WinBUGS object to an rv object}
\alias{as.rv.bugs}
\usage{
  as.rv.bugs(x)
}
\description{
  \code{as.rv.bugs} coerces an \code{R2WinBUGS} object to an \code{rv} object.
}
\arguments{
  \item{x}{a bugs (R2WinBUGS) object}
}
\details{
 
}
\value{
  A random vector with the names attribute set 
  to the names of the parameters.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\seealso{
  \code{\link{split.rv}} to split the named random vector so you can
  ``attach" it with \code{\link{rvattach}}.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}

\keyword{classes}

