% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitbyname.R
\name{splitbyname}
\alias{splitbyname}
\title{Split a vector based on the names of the components}
\usage{
splitbyname(x)
}
\arguments{
\item{x}{a vector or a list with the name attributes set}
}
\value{
A list of arrays and vectors.  Missing entries in the arrays and
vectors are filled in with \code{NA}s.
}
\description{
\code{splitbyname} is a utility function that splits the given vector based
on the names of the components and returns a named list of arrays and
vectors.
}
\details{
The names are supposed to be of the format \code{name[index]}, for example
\code{alpha[1,1]}, \code{beta[1]}, etc.

A name without brackets is equivalent to a name with \verb{[1]}.

The dimension attribute will not be set in case of vectors.
}
\examples{

  x <- structure(c(1,3), names=c("x[1,1]", "x[3,3]"))
  splitbyname(x) # yields a list containing a 3x3 matrix

}
\author{
Jouni Kerman \email{jouni@kerman.com}
}
\keyword{manip}
