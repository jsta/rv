\name{rvfactor}
\title{Categorical Random Variables (Random Factors)}
\alias{rvfactor}
\alias{rvfactor.rv}
\alias{is.rvfactor}
\alias{as.rvfactor}
\alias{as.rv.rvfactor}
\usage{
  rvfactor(x, \dots)
  \method{rvfactor}{rv}(x, levels=NULL, \dots)
  is.rvfactor(x)
  \method{as.rv}{rvfactor}(x)
  as.rvfactor(x, \dots)
}
\description{
  Creates or tests for objects of type ``\code{rvfactor}".
}
\arguments{
  \item{x}{object to be coerced or tested.}
  \item{levels}{factor levels (labels for the levels)}
  \item{\dots}{other arguments}
}
\details{
  Internally random factors are integer-valued just like regular factors in R.
}
\value{
  \code{rvfactor}: an \code{rvfactor} object.

  \code{is.rvfactor}: \code{TRUE} or \code{FALSE}.

  \code{as.rv.rvfactor}: an \code{rv} object.

  \code{as.rvfactor.rv}: an \code{rvfactor} object.
}
\note{
}
\seealso{
  
}
\references{Kerman, Jouni and Gelman, Andrew. 
Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
\url{http://www.stat.columbia.edu/~kerman/Research/postsim.pdf}
}

\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}

\examples{
  # Probabilities of each integer of trunc(Z) where Z ~ N(0,1) ?
  x <- rvnorm(1)
  rvfactor(trunc(x))
  rvfactor(x>0)
  rvfactor(rvpois(1, lambda=0.5))
}
\keyword{classes}

