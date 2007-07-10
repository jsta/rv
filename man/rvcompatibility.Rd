\name{rvcompatibility}
\title{Set Compatibility Level of the rv Package}
\alias{rvcompatibility}
\usage{
  rvcompatibility(level=NULL)
}
\description{
  Sets the compatibility level of the rv package.
}
\arguments{
  \item{level}{0, 1, or NULL. 0 or 1 sets the compatibility level; NULL causes the function to return the current level.} 
}
\details{
  'rv' implements (loosely) a \emph{replacement} for the standard numeric and logical classes.
  A numeric or logical vector is just a special case of the more general random vector data type (same would apply for the complex data type but this isn't implemented).
  To implement the package right, we need certain functions to work completely with
  random vectors and basically 'fool' R into believing that the 'rv' objects are 
  just as good as the standard R data types. This is not possible, unfortunately, 
  without replacing some of the "primitive" non-generic functions. 
  Everything would be fine if everything were generic and we could re-define 
  what 'numeric' means for example. This isn't possible (for obvious efficiency 
  and compatibility reasons I think) so we need to do some surgery, replacing
  standard functions with special ones that check whether a given argument is 
  an rv object or not, and then choose the right function to work on the arguments
  (a primitive kind of a method dispatch mechanism). 
  Unfortunately even this will not make rv fully compatible. 
  Especially, the function \code{"[<-"}, even when replaced, will not work with 
  all programs. Thus we need a system to undo a function replacement. 
  If you experience a weird error (esp. working with assignment to a factor vector)
  check \code{rvcompatibility(1)} and re-run your program.

  The only use of this function - currently - is \code{rvcompatibility(1)} 
  which restores the function \code{"[<-"} in the base package.
  Inserting numeric/logical objects into an rv objects still works, 
  but inserting an rv object into a plain numeric/logical vector will not work.

  NOTE. When the package rv is detached (you can do it simply with \code{detachrv()}),
  all functions that were replaced upon attaching the package are immediately restored.
  
  \item{level 0}{Default: In particular, makes sure that assignment of an rv object into a component of numeric object is allowed (returning an rv object)}
  \item{level 1}{More compatible with existing programs. In particular, assignment of an rv object into a component of numeric object is not possible, which is a default behaviour of R since R does not allow }
}
\value{
  If \code{level} is \code{NULL}, the current compatibility level.

  If \code{level} is 0 or 1, returns the \emph{previous} compatibility level.
}
\examples{
 \dontrun{
  rvcompatibility(0)
  x <- 1:9
  # Will work since rv checks that the value to be injected is an rv:
  x[1] <- rvnorm(1)
  print(x)
  \dontrun{
    rvcompatibility(1)
    x <- 1:9
    # Will fail since R thinks this is assigning a list into a vector component:
    x[1] <- rvnorm(1)
    print(class(x)) # "list"
  }
  # However, this works:
  x <- as.rv(1:9)
  rvcompatibility(1)
  x[1] <- rvnorm(1)
  rvcompatibility(0) # restore the fully rv-enabled computing environment...
 }
}
\keyword{classes}
