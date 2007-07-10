\name{logical.rv}
\title{Logical Random vectors}
\alias{logical.rv}
\alias{as.logical.rv}
\alias{is.logical.rv}
\usage{
  \method{is.logical}{rv}(x)
  \method{as.logical}{rv}(x, \dots)
}
\description{
  \code{is.logical(x)} returns \code{TRUE} for each component that 
  contains an logical-valued vector, otherwise \code{FALSE} is given.
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{Further arguments passed on}
}
\details{
  \code{is.logical(x)} returns \code{TRUE} if and only if \emph{each} component 
  of \code{x} is logical-valued (i.e. \code{TRUE}/\code{FALSE}).

   \code{as.logical.rv} coerces an rv object into logical-valued one.
   In effect, the function  \code{as.logical} is applied to all simulations.

}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\seealso{
  \code{\link{as.logical.rv}},  \code{\link{is.logical.rv}}
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvbern(prob=0.5)   # some 0/1 valued random variable
  is.logical(x)           # FALSE, because by default x is 'double'
  x <- as.logical(x)      # coerce to logical; all zeros become FALSE, ones become TRUE
  is.logical(x)           # TRUE
}

\keyword{classes}

