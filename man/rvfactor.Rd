\name{rvfactor}
\title{Categorical Random Variables (Random Factors)}
\alias{rvfactor}
\alias{rvfactor.rv}
\alias{is.rvfactor}
\alias{as.rvfactor}
\alias{as.rv.rvfactor}
\alias{print.rvfactor}
\usage{
  rvfactor(x, \dots)
  \method{rvfactor}{rv}(x, levels=NULL, \dots)
  is.rvfactor(x)
  \method{as.rv}{rvfactor}(x, \dots)
  as.rvfactor(x, \dots)
  \method{print}{rvfactor}(x, all.levels=FALSE, \dots)
}
\description{
  Creates or tests for objects of type ``\code{rvfactor}".
}
\arguments{
  \item{x}{object to be coerced or tested.}
  \item{levels}{factor levels (labels for the levels)}
  \item{all.levels}{logical; whether to print all levels or not (see below for details)}
  \item{\dots}{other arguments}
}
\details{
  Internally random factors are integer-valued just like regular factors in R.

  The number of levels to print when \code{all.levels==FALSE} 
  can be set by \code{rvpar(max.levels=...)}. By default this is set to 10.
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
  # Probabilities of each integer of trunc(Z) where Z ~ N(0,1) ?
  x <- rvnorm(1)
  rvfactor(trunc(x))
  rvfactor(x>0)
  rvfactor(rvpois(1, lambda=0.5))
}
\keyword{classes}

