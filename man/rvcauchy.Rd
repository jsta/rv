\name{rvcauchy}
\title{Generate Random Variables from a Cauchy Sampling Model}
\alias{rvcauchy}
\usage{
rvcauchy(n=1, location=0, scale=1)
}
\description{
  Random vector generation for the Cauchy distribution.
}
\arguments{
  \item{n}{integer: number of variables to generate}
  \item{location}{location parameter (may be random)}
  \item{scale}{scale parameter (may be random)}
}
\details{
  For details on the Cauchy distribution, see \link{Cauchy}.
  See also \code{\link{rvt}}; Cauchy is a special 
  case of the t-distribution with 1 degree of freedom,
  and therefore \code{rvcauchy(n,location,scale)} is 
  equivalent to \code{rvt(n, mu, scale, df=1)}.
}
\value{
  A random vector (rv object) of length \code{n}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}

\keyword{classes}

