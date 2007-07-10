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
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{classes}

%%%
