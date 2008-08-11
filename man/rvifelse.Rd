\name{rvifelse}
\title{Conditional Random Element Selection}
\alias{rvifelse}
\usage{
  rvifelse (test, yes, no)
}
\description{
  \code{rvifelse} is the \code{rv}-compatible version of the function \code{ifelse}.
}
\arguments{
  \item{test}{an object which can be coerced to logical mode.}
  \item{yes}{return values for true elements of \code{test}}
  \item{no}{return joint simulations and not simulations from each component separately}
}
\details{
  \code{rvifelse} returns a \emph{random} value with the same shape as
  \code{test} which is filled with random or constant elements selected from either
  \code{yes} or \code{no},
  depending on whether the random draw in an
  element of \code{test} is \code{TRUE} or \code{FALSE}.

  
}
\value{
  A \emph{numeric} array of dimensions \code{size} times \code{length(x)}.
}
\seealso{
  \code{\link{ifelse}}.
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  
}

\keyword{classes}

