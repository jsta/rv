\name{rvmixed}
\title{Categorical Random Variables (Random Factors)}
\alias{rvmixed}
\alias{rvmixed.rv}
\alias{is.rvmixed}
\alias{is.rvmixed.rv}
\alias{as.rv.rvmixed}
\usage{
  rvmixed(x, pointmass=NULL)
  is.rvmixed(x)
  \method{is.rvmixed}{rv}(x)
  \method{as.rv}{rvmixed}(x)
}
\description{
  Creates or tests for objects of type ``\code{rvmixed}".
}
\arguments{
  \item{x}{object to be coerced or tested.}
  \item{pointmass}{vector of points with point mass probability}
  \item{\dots}{other arguments}
}
\details{
  Internally random factors are integer-valued just like regular factors in R.
}
\value{
  \code{rvmixed}: an \code{rvmixed} object.

  \code{is.rvmixed}: \code{TRUE} or \code{FALSE}.

  \code{as.rv.rvmixed}: an \code{rv} object.

  \code{as.rvmixed.rv}: an \code{rvmixed} object.
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
  #
}
\keyword{classes}

