\name{rvattach}
\title{Attach components of rv objects to search path}
\alias{rvattach}
\usage{
rvattach(what, name='rvattach', overwrite=TRUE, impute=FALSE, \dots)
}
\description{
  
  \code{rvattach} 
}
\arguments{
  \item{what}{an rv object to be split up into a list of sub-vectors and attached}
  \item{name}{}
  \item{overwrite}{If \code{TRUE}, objects with identical names in the Workspace (.GlobalEnv) that are masking objects in the database to be attached will be deleted.}
  \item{impute}{If \code{TRUE}, the components of each sub-vector are imputed into the object with identical name in the Workspace and saved into the list. The object in the Workspace is not affected.}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  \code{rvattach} takes an rv object, splits it up into a list
  of sub-vectors 
  (using \code{split.rv}) by their name attributes and attaches the list
  to the search path. 
  For example, if the object to be processed has names \code{x[1]} and \code{y[3]},
  this will be split into a list of two sub-vectors, \code{x} and \code{y}, 
  and then attached.
  If the option impute is \code{TRUE}, the corresponding variables
  with the same names are merged with the ones that are found in .GlobalEnv. 
  If the option overwrite is \code{TRUE}, the objects with identical
  names in the Workspace are deleted.
}
\references{Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{classes}

%%%%
