\name{splitbyname}
\title{Split a vector based on the names of the components}
\alias{splitbyname}
\usage{
  splitbyname(x)
}
\description{
  \code{splitbyname} splits the given vector based on the names
  of the components and returns a named list of
  arrays and vectors.
}
\arguments{
  \item{x}{a vector with the name attributes set}
}
\details{
  The names are supposed to be of the format 'name[indexes]',
  for example 'alpha[1,1]', 'beta[1]', etc.

  A name without brackets is equivalent to a name with '[1]'.

  Dimension will not be set in case of vectors.

  Note. Will not work for rv objects if \code{rvcompatibility()==1}.
}
\value{
  A list of arrays and vectors.
  Missing entries in the arrays and vectors are filled in with \code{NA}s.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- structure(rvnorm(2), names=c("x[1,1]", "x[3,3]"))
  splitbyname(x) # yields a list containing a 3*3 random matrix
}

\keyword{manip}

%%%%
