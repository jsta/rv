\name{rvAdapted}
\title{Functions adapted to accept random variable (rv) objects}
\alias{rvAdapted}
\alias{c}
\alias{\%*\%}
\alias{is.vector}
\alias{max}
\alias{pmax}
\alias{min}
\alias{pmin}
\alias{median}
\alias{plot}
\usage{
 c(\dots, recursive=FALSE)
 x \%*\% y
 is.vector(x, mode="any")
 max(\dots, na.rm = FALSE)
 pmax(\dots, na.rm = FALSE)
 min(\dots, na.rm = FALSE)
 pmin(\dots, na.rm = FALSE) 
 median(x, na.rm = FALSE) 
}
\description{
Functions that have been adapted to work with rv objects.
}
\details{
  These functions have been adapted to work with rv objects
  (since no \emph{real} object-oriented generic version is available).
  They check whether an rv object has been passed as an argument;
  if false, the original function will be called. 

  To restore the original functions, detach the package; 
  use \code{detachrv()}.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\seealso{
  \code{\link[rv]{cov}} for details about the \code{cov} and \code{var}
  functions adapted for random variable objects.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  x <- rvnorm(10) # A random vector.
  median(x)       # A random variable.
  max(x)          # Another random variable.  
}

\keyword{internal}

%  The original help pages are available directly via these links:
%  \code{\link[base]{c}},
%  \code{\link[base]{\%*\%}},
%  \code{\link[base]{is.vector}},
%  \code{\link[base]{max}},
%  \code{\link[base]{min}},
%  \code{\link[base]{pmin}},
%  \code{\link[base]{pmax}},
%  \code{\link[base]{median}},
%  \code{\link[graphics]{plot}}.
