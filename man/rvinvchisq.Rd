\name{rvinvchisq}
\title{Generate Random Variables from a Inverse-Chi-Square Sampling Model}
\alias{rvinvchisq}
\usage{
  rvinvchisq(n=1, df, scale=1)
}
\description{
  \code{rvinvchisq}
}
\arguments{
  \item{n}{integer: number of variables to generate}
  \item{df}{degrees of freedom (may be random)}
  \item{scale}{scale parameter (may be random)}
}
\details{
  \code{rvinvchisq} 
}
\value{
  A random vector (rv object).
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
  rvinvchisq(df=3, scale=2)
}

\keyword{classes}

