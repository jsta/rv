\name{rvnsims}
\title{Set or Get the Default Number of Simulations to Generate for rv Objects}
\alias{rvnsims}
\usage{
rvnsims(nsims)
}
\description{
  Sets or retrieves the default number of simulations.
}
\arguments{
  \item{nsims}{Number of simulations.}
}
\details{
  %...
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}

\examples{
  rvnsims(10000)
  print(rvnsims()) # prints 10000
}


\keyword{classes} %??
