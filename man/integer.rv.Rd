\name{integer.rv}
\title{Integer-valued Random Vectors}
\alias{integer.rv}
\alias{is.integer.rv}
\alias{as.integer.rv}
\usage{
  \method{is.integer}{rv}(x)
  \method{as.integer}{rv}(x, \dots)
}
\description{
  Tests or coerces (numeric or logical-valued) 
  random vector objects into integer-valued ones.
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{other arguments}
}
\details{
  \code{is.integer(x)} returns \code{TRUE} if and only if \emph{each} component 
  of \code{x} is integer-valued.

   \code{as.integer.rv} coerces an rv object into integer-valued one.
   In effect, the function  \code{as.integer} is applied to all simulations.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\seealso{
  \code{\link{logical.rv}}
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvpois(lambda=1)   # some integer-valued random variable
  is.integer(x)           # FALSE, because by default x is 'double'!
  x <- as.integer(x)      # coerce to integer
  is.integer(x)           # TRUE
}

\keyword{classes}

