% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvsimapply.R
\name{rvsimapply}
\alias{rvsimapply}
\title{Apply a Function to Columns of the Matrix of Simulation of a Random Vector}
\usage{
rvsimapply(x, FUN, ...)
}
\arguments{
\item{x}{an object}

\item{FUN}{an R function object}

\item{\dots}{further arguments passed to the function \code{FUN}}
}
\value{
A numeric vector or array.
}
\description{
\code{rvsimapply}
}
\details{
\code{rvsimapply} applies a given function to the \emph{rows} of the
simulation matrix of the given random vector.

If the function is to be applied to \emph{rows} of the simulation matrix,
use \code{\link{simapply}} or \code{\link{rvmapply}} instead.

Usually used in functions that implement an 'rv'-compatible routine.
}
\examples{

  #

}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects. Statistics and
Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\author{
Jouni Kerman \email{jouni@kerman.com}
}
\keyword{classes}
