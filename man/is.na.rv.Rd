\name{is.na.rv}
\title{Missing Data Indicators}
\alias{is.na.rv}
\alias{rv.all.na}
\alias{rv.any.na}
\usage{
  \method{is.na}{rv}(x)
  rv.all.na(x)
  rv.any.na(x)
}
\description{
  \code{is.na.rv} returns the distribution (random variable)
  of the indicator function of missing data.
  \code{rv.all.na} returns \code{TRUE} if all components
  of the argument vector are completely missing.
  \code{rv.any.na} returns \code{TRUE} if any component
  of the argument vector has missing values.
}
\arguments{
  \item{x}{an rv object}
}
\details{
  Internally, \code{is.na.rv} applies the function \code{is.na}
  to each simulation of each component of the argument vector.
}
\value{
 \code{is.na.rv} returns a ``Bernoulli'' random vector 
  of the same length and dimension as those of \code{x}.

 \code{rv.all.na} and \code{rv.any.na} return \code{TRUE} or \code{FALSE} (single value).
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
  x <- trunc(rvnorm(1))
  y <- !(x==0 & NA) # TRUE if x!=0
  x <- y*x
  is.na(x)     # 69%: Pr(-1<Z<1)
  is.logical(is.na(x)) # TRUE
  rv.any.na(x) # TRUE
  rv.all.na(x) # FALSE
}

\keyword{classes}

%%% maybe make up an example where we have random-length Poisson variable?
