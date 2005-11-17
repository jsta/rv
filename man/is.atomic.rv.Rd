\name{is.atomic.rv}
\title{Random vectors are atomic}
\alias{is.atomic.rv}
\usage{
  \method{is.atomic}{rv}(x)
}
\description{
  \code{is.atomic(x)} returns \code{TRUE} for rv objects.
  This is necessary to make R believe that rvs are as good as other vectors.
}
\arguments{
  \item{x}{an object}
}
\details{
  \code{is.atomic.rv} 
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  is.atomic(rvnorm(1)) # TRUE, believe it or not.
}

\keyword{classes}

