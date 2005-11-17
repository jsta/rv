\name{detachrv}
\title{Detach the rv package}
\alias{detachrv}
\usage{
detachrv()
}
\description{
  \code{detachrv} detaches the rv package and restores the original
  functions in base and stats packages.
  Currently equivalent to \code{detach("package:rv")}.
}
%\arguments{
%  
%}
\details{
  \code{detachrv} 
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  \dontrun{
  library(rv)
  detachrv()
  }
}

\keyword{classes}


