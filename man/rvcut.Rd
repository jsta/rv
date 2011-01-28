\name{rvcut}
\title{Convert Numeric to Random Factor}
\alias{rvcut}
\alias{rvcut.rv}
\usage{
  rvcut(x, \dots)
  \method{rvcut}{rv}(x, \dots)
}
\description{
  Convert implements the `cut' function using random variables.
}
\arguments{
  \item{x}{a plain or a random vector which is to be converted to a factor by cutting.}
  \item{\dots}{arguments passed to the function \code{\link{cut}}.}
}
%\details{
%    
%}
\value{
  A random factor.
}
\seealso{
  \code{\link{rvfactor}}, \code{\link{cut}}.
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
  rvcut(rvnorm(1), breaks=c(-Inf,-2,-1,0,1,2,Inf))
}
\keyword{classes}

