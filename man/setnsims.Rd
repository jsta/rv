\name{setnsims}
\title{Set or Get the Default Number of Simulations to Generate for rv Objects}
\alias{setnsims}
\alias{getnsims}
\usage{
  setnsims(nsims)
  getnsims()
}
\description{
  Sets or retrieves the default number of simulations.
}
\arguments{
  \item{nsims}{Number of simulations.}
}
\details{
  \code{getnsims()} is an alias for the function call \code{setnsims()}.
}
\value{
  If \code{nsims} is specified, returns the \emph{old} value, 
  else returns the current value.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}

\examples{
  setnsims(1000) # set n.sims to 1000
  s <- setnsims(10000) # s is now 1000
  print(getnsims()) # prints 10000
  setnsims(s) # restore the number of simulations back to 1000
}


\keyword{classes} %??
