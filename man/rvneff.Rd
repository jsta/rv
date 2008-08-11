\name{rvneff}
\title{Number of Effective Draws in Each Component of a Random Variable}
\alias{rvneff}
\usage{
  rvneff(x)
}
\description{
  Retrieves the number of effective draws in each component of the argument.
}
\arguments{
  \item{x}{an rv object}
}
\details{
  The number of effective draws is supposed to be saved by the simulation generating program (e.g. WinBUGS via R2WinBUGS).
}
\value{
  A numeric object of the same length as the argument \code{x}.
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
  #
}

\keyword{classes}

%%%
