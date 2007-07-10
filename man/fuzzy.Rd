\name{fuzzy}
\title{Fuzziness}
\alias{fuzzy}
\alias{is.fuzzy}
\alias{is.fuzzy.rv}
\usage{
  is.fuzzy(x)
  \method{is.fuzzy}{rv}(x)
}
\description{
  Tests whether an object is ``fuzzy",
  i.e.
  a logical random scalar 
  that has probability strictly between zero and one (not strictly true nor strictly false).
}
\arguments{
  \item{x}{an object, random or constant}
}
\details{
    
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- as.logical(rvbern(1,0.4)) # a logical random variable
  is.fuzzy(x) # TRUE, since x is logical and not constant
  is.fuzzy(x<2) # FALSE, since x is less than 2 with probability one
  is.fuzzy(rvnorm(1)) # FALSE, since it's not a probability
  is.fuzzy(TRUE) # FALSE, since TRUE is strictly TRUE
  is.fuzzy(1) # FALSE, since 1 is not a logical variable
}

\keyword{classes}

