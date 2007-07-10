\name{sort.rv}
\title{Distribution of Order Statistics of a Random Vector}
\alias{sort.rv}
\usage{
  \method{sort}{rv}(x, \dots)
}
\description{
  \code{sort.rv} computes the distribution of the order statistics
  of a random vector.
}
\arguments{
  \item{x}{a random vector}
  \item{\dots}{further arguments passed to \code{sort.rv}}
}
\details{
 The result is the \emph{distribution} of the order statistic of the given vector \code{x}:
 that is, the \code{sort} function is applied to each \emph{row}
 of the matrix of simulations of \code{x} (\code{sims(x)}) and returned then
 in random vector form. 
 
 See  \code{\link{sort}} for further details of the function \code{sort}.
}\value{
  An rv object of the same length as \code{x}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\seealso{
  \code{\link{sort}}
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{manip}

