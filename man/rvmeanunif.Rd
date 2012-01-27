\name{rvmeanunif}
\alias{rvmeanunif}
\alias{rvtriang}
\title{
 The distribution of the mean of uniform random variables
}
\description{
  The distribution of the mean of uniform random variables
  with each of them in the interval \code{(-1, 1)}.
}
\usage{
  rvtriang(n = 1, mode, scale)
  rvmeanunif(n=1, mode, scale, df)
}
\arguments{
  \item{n}{Length of the vector to output}
  \item{mode}{Mode (center) of the distribution}
  \item{scale}{Scale (half-width) of the distribution around the mode}
  \item{df}{``degrees of freedom'': number of independent components to
    average}
}
\details{
  Assuming that all inputs are constants, 
  each generated variable has a mode (center) at \code{mode},
  constrained between \code{(-scale, scale)}. 

  The shape becomes more and more bell-shaped (Normal) as the number of
  the independent variables in the sum (mean) increases.

  The case of \code{df=2} (mean of two variables) is the special case of
  the symmetric triangular distribution in the range 
}
\value{
  A random vector of length \code{n}.
}
\author{
  J Kerman
}
\examples{
  x <- rvtriang(1)
  y <- rvmeanunif(df=2) ## same distribution as that of x
}
\keyword{ dist }


