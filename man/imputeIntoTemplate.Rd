\name{imputeIntoTemplate}
\title{Force a vector into specified shape}
\alias{imputeIntoTemplate}
\usage{
  imputeIntoTemplate(x, template)
}
\description{
  \code{imputeintoTemplate} imputes 'x' into the array shape given by 'template',
  which is a logical array of the desired dimension
  (usually the dimnames are also important!)
}
\details{
  The template is a list with names matching the names of the components of '\code{x}'.
  The components of the list are vectors or arrays of logical values (and NAs) 
  with (typically) dimnames set.

  The values in the components of the template vector/array
  determine how the imputation is to be done.

  A 'TRUE' value indicates a component that is replaced
  by a value in x.

  Value 'FALSE' indicates a component that is NOT imputed,
  but skipped and later imputed an 'NA'.

  A component with 'NA' in the template \emph{is imputed} but
  it is later masked with an 'NA'.

  Thus only components with the value 'TRUE' will appear non-NA
  in the return value.
}
\arguments{
  \item{x}{a vector with the name attributes set}
  \item{template}{a named list of templates}
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  #
}

\keyword{manip}

%%%%
