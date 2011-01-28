\name{as.logical.rv}
\title{Logical Random vectors}
\alias{as.logical.rv}
\usage{
  \method{as.logical}{rv}(x, \dots)
}
\description{
  Coerces a random variable to a logical-valued one (Bernoulli r.v.)
}
\arguments{
  \item{x}{an rv object}
  \item{\dots}{Further arguments passed on}
}
\details{
   In effect, the function  \code{as.logical} is applied to all simulations.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
%\seealso{
%}
\note{
  \code{is.logical(x)} returns \code{TRUE} if and only if \emph{each} component 
  of \code{x} is logical-valued (i.e. \code{TRUE}/\code{FALSE}).
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  x <- rvbern(prob=0.5)   # some 0/1 valued random variable
  print(x)
  is.logical(x)           # FALSE, because by default x is 'double'
  x <- as.logical(x)      # coerce to logical; all zeros become FALSE, ones become TRUE
  is.logical(x)           # TRUE
  print(x)                # Shows the expectations and not the quantiles
}

\keyword{classes}

