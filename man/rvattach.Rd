\name{rvattach}
\title{Attach components of rv objects to search path}
\alias{rvattach}
\alias{rattach}
\usage{
  rattach(what, name=deparse(substitute(what)), overwrite=TRUE, impute=FALSE, \dots)
  rvattach(what, name=deparse(substitute(what)), overwrite=TRUE, impute=FALSE, \dots)
}
\description{
  \code{rattach} works like \code{attach} but removes 
  objects from \code{.GlobalEnv} that have conflicting names;
  optionally can also \emph{impute} vector components \emph{by name} into the 
  object that is in \code{.GlobalEnv} (copying the imputed object into the
  attached list.
  \code{rvattach} first splits a named vector into a list, and then 
  uses \code{rattach} to attach them to the search path.
}
\arguments{
  \item{what}{an rv object to be split up into a list of sub-vectors and attached}
  \item{name}{name to use for the attached database. All databases with the exact same name in the search path are removed first.}
  \item{overwrite}{If \code{TRUE}, objects with identical names in the workspace (\code{.GlobalEnv}) that are masking objects in the database to be attached will be deleted.}
  \item{impute}{If \code{TRUE}, the components of each sub-vector are imputed into the object with identical name in the workspace and saved into the list. The object in the workspace is not affected (unless \code{overwrite=TRUE}).}
  \item{\dots}{further arguments passed to \code{attach}}
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
  with the same names are merged with the ones that are found in \code{.GlobalEnv},
  and copied into the database.

  If the option overwrite is \code{TRUE}, the objects with identical
  names in the Workspace are deleted.

  All databases with the exact same name in the search path are removed first.
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
  # Trying rattach:
  x <- rvnorm(9)
  names(x) <- paste("theta[", 1+seq_along(x), "]", sep="")
  L <- list(x=x)
  x <- 1:10
  rattach(L, impute=TRUE) # attaches a list called 'L' containing an object 'x'
  print(x)
  detach("L")
  exists("x") # FALSE
  # Trying rvattach:
  rvattach(L$x) # attaches a list 'rvattach' containing an object 'theta'
  print(theta)
  rvattach() # Detaches the list 'rvattach'
  exists("theta") # FALSE
}

\keyword{classes}

%%%%
