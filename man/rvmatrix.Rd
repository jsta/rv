\name{rvmatrix}
\title{Matrices and Arrays of Random Vectors}
\alias{rvmatrix}
\alias{rvarray}
\usage{
 rvmatrix(\dots)
 rvarray(\dots)
}
\description{
  Arrange a given random vector into a matrix or array form.
}
\arguments{
  \item{\dots}{arguments passed to \code{matrix} or to \code{array}}
}
\details{
  These are 'rv' compatible versions of the functions 
  \code{\link{rvmatrix}} and  \code{\link{rvarray}}.
}
\seealso{
  \code{\link{mlplot}}
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  n.rows <- 3; n.cols <- 4; n <- (n.rows*n.cols)
  mu.true <- rnorm(1:n.rows, mean=1:n.rows, sd=1)
  theta <- rvmatrix(rvnorm(n=n.cols, mean=mu.true, sd=0.5), nrow=n.rows)
  col.labels <- paste("Time", 1:n.cols, sep=":")
  row.labels <- paste("Unit", 1:n.rows, sep=":")
  dimnames(theta) <- list(row.labels, col.labels)
  print(theta)
  print(E(theta))  
}

\keyword{classes}

