\name{splitAndImpute}
\title{Split a vector based on the names of the components}
\alias{splitAndImpute}
\usage{
  splitAndImpute(x, templates)
}
\description{
  \code{splitAndImpute} splits the given vector based on the names
  of the components, and imputes the variables into given templates.
}
\arguments{
  \item{x}{a vector with the name attributes set}
  \item{templates}{a named list of templates (arrays)}
}
\details{
  A convenience function: combination of
  \code{splitbyname} and \code{imputeIntoTemplate}.

  See \code{\link{splitbyname}} and \code{\link{imputeIntoTemplate}} for details.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{manip}

%%%%
