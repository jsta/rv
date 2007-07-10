\name{numeric.rv}
\title{Numeric Random Vectors}
\alias{numeric.rv}
\alias{as.numeric.rv}
\alias{is.numeric.rv}
\alias{is.numeric.rvfactor}
\alias{is.numeric.rvmixed}
\usage{
  \method{is.numeric}{rv}(x)
  \method{as.numeric}{rv}(x, \dots)
  \method{is.numeric}{rvfactor}(x)
  \method{is.numeric}{rvmixed}(x)
}
\description{
  \code{is.numeric(x)} returns \code{TRUE} if and only if \emph{all} components
  of an rv object are numeric.
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{Further arguments passed on}
}
\details{
  \code{is.numeric(x)} returns \code{TRUE} if and only if \emph{each} component 
  of \code{x} is numeric-valued (i.e. \code{TRUE}/\code{FALSE}).

   \code{as.numeric.rv} coerces an rv object into numeric-valued one.
   In effect, the function  \code{as.numeric} is applied to all simulations.

  Random factors are not numeric; mixtures are. 
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
  x <- as.logical(rvbern(1,0.5)) # Bernoulli rv
  is.numeric(x)           # FALSE
  x <- as.numeric(x)      # coerce to numeric; all TRUEs become ones, FALSEs zeros
  is.numeric(x)           # TRUE
}

\keyword{classes}

