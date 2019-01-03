#' Extract or Replace Parts of a Random Vector
#' 
#' @name Extract-rv
#' 
#' @description Bracket slice and assignment methods adapted for random vectors and arrays.
#' The assignment function \code{impute<-} is compatible with both non-rv and
#' rv objects (rv, rvsummary, and rvfactor objects). To write universal code
#' that works both atomic and rv objects, use \code{impute(x, ...) <- value}
#' instead of \code{x[...] <- value}.
#' 
#' 
#' NOTE. \code{x} will NOT be automatically coerced into an rv object.
#' 
#' \code{value} may be an rv object or a regular numeric object.
#' 
#' Extracting rv objects works the same way as extracting components of a
#' numerical vector or array.  The return value is always an object of class
#' 'rv'.  Type ?Extract for details.
#' 
#' Note: the index arguments (\code{i}, \code{j}, etc.)  \emph{must} be
#' constants, but this may change in the future.
#' 
#' %Note: the index arguments (\code{i}, \code{j}, etc.) may be %themselves
#' random variables, however they will be coerced %into \emph{integers}, as one
#' would expect.
#' 
#' @aliases [.rv [.rvfactor [.rvsummary [<-.rv [<-.rv [<-.rvsummary impute<-
#' @param x object from which to extract element(s) or in which to replace
#' element(s).
#' @param \dots indices specifying elements to extract or replace.
#' @param value typically an array-like R object of a similar class as
#' \code{x}.
#' @param drop For matrices and arrays.  If \code{TRUE} the result is coerced
#' to the lowest possible dimension (see the examples).  This only works for
#' extracting elements, not for the replacement.
#' @return A random variable (an rv object).
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(1)
#'   y <- (1:5)
#'   \dontrun{
#'      y[2] <- x ## Will not work
#'   }
#'   impute(y, 2) <- x
#' 
NULL


#' @export 
"[.rv" <- function (x, ..., drop = TRUE)
{
  cx <- class(x)
  X <- NextMethod()
  class(X) <- cx
  return(X)
}

#' @export 
"[.rvsummary" <- function (x, ..., drop = TRUE)
{
  q <- attr(x, "quantiles")
  cx <- class(x)
  x <- NextMethod()
  class(x) <- cx
  attr(x, "quantiles") <- q
  return(x)
}

#' @export 
"[<-.rvsummary" <- function (x, ..., value = NULL)
{
  cx <- class(x)
  q <- attr(x, "quantiles")
  value <- as.rvsummary(value, quantiles=q)
  X <- .Primitive("[<-")(unclass(x), ..., value=value)
  class(X) <- cx
  return(X)
}

#' @export
"[<-.rv" <- function (x, ..., value = NULL)
{
  cx <- class(x)
  value <- as.rvobj(value)
  X <- .Primitive("[<-")(unclass(x), ..., value=value)
  class(X) <- cx
  return(X)
}

#' @inheritParams Extract-rv
#' @export
"impute<-" <- function(x, ..., value) {
  if (! is.rvobj(x) && ! is.rvobj(value)) {
    x[...] <- value
  } else if (is.rvsummary(x) || is.rvsummary(value)) {
    x <- as.rvsummary(x)
    value <- as.rvsummary(value)
    x[...] <- value
  } else if (is.rv(x) || is.rv(value)) {
    x <- as.rv(x)
    value <- as.rv(value)
    x[...] <- value
  }
  return(x)
}
