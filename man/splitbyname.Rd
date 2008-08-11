\name{splitbyname}
\title{Split a vector based on the names of the components}
\alias{splitbyname}
\usage{
  splitbyname(x)
}
\description{
  \code{splitbyname} is a utility function that 
  splits the given vector based on the names
  of the components and returns a named list of
  arrays and vectors.
}
\arguments{
  \item{x}{a vector or a list with the name attributes set}
}
\details{
  The names are supposed to be of the format 'name[index]',
  for example 'alpha[1,1]', 'beta[1]', etc.

  A name without brackets is equivalent to a name with '[1]'.

  The dimension attribute will not be set in case of vectors.

}
\value{
  A list of arrays and vectors.
  Missing entries in the arrays and vectors are filled in with \code{NA}s.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  x <- structure(c(1,3), names=c("x[1,1]", "x[3,3]"))
  splitbyname(x) # yields a list containing a 3x3 matrix
}

\keyword{manip}

%%%%
