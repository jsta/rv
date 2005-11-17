\name{Summary.rv}
\title{Summary of a Random Vector}
\alias{Summary.rv}
\usage{
Summary.rv(\dots, na.rm=FALSE)
}
\description{
  \code{Summary.rv}
}
\arguments{
  \item{\dots}{rv object}
  \item{na.rm}{logical, remove NAs?}
}
\details{
  
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvnorm(6)
  dim(x) <- c(2,3)
  summary(x)
}

\keyword{classes}

