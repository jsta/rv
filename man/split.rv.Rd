\name{split.rv}
\title{Split a Random Vector into a List by Name}
\alias{split.rv}
\usage{
  split.rv(x, f, drop=FALSE, ...)
}
\description{
 Split a random vector into a list by name.
}
\arguments{
  \item{x}{a random vector}
  \item{f}{a ``factor" that specifies the grouping, by default the variable names derived from \code{names{x}}.}
  \item{drop}{logical indicating if levels that do not occur should be dropped (if f is a factor or a list).}
  \item{\dots}{further potential arguments passed to methods.}
}
\details{
  If \code{f} is omitted, splitting is done by \code{names(x)}, but
  with the possible brackets (as in ``\code{alpha[10]}") omitted.
}
\seealso{
  \code{\link{split}}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvnorm(3)
  names(x) <- c("alpha[1]", "alpha[2]", "theta")
  split(x)
}

\keyword{classes}

