#' Generate a Random Vector from a Bernoulli Sampling Model
#' 
#' \code{rvmatch} returns a (random) vector of the positions of (first) matches
#' of its first argument in its second.
#' 
#' \code{\%*in*\%} is a binary operator (analogous in its operation to
#' \code{\%in\%}) which returns a logical (random) vector indicating if there is
#' a match or not for its left operand.
#' 
#' ...
#' 
#' @aliases rvmatch 
#' @param x random vector, regular atomic vector, or \code{NULL}: the values to
#' be matched.
#' @param table random vector, regular atomic vector, or \code{NULL}: the
#' values to be matched against.
#' @param nomatch the value to be returned in the case when no match is found.
#' Note that the value is coerced to \code{integer}.
#' @param incomparables a vector of values that cannot be matched.  Any value
#' in x matching a value in this vector is assigned the nomatch value.  For
#' historical reasons, FALSE is equivalent to \code{NULL}
#' @return A random vector (an rv object) of the same length as \code{x}.
#' 
#' \code{rvmatch} returns an integer-valued vector.
#' 
#' \code{\%*in*\%} returns a logical-valued vector.
#' 
#' Both functions are compatible with regular atomic vectors.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvempirical(5, 1:10)
#'   z <- rvmatch(x, table=1:3, nomatch=0L)
#' 
#' @export rvmatch
rvmatch <- function (x, table, nomatch=NA_integer_, incomparables=NULL) {
  rvmapply(base::match, x=x, table=table, nomatch=nomatch, MoreArgs=list(incomparables=incomparables))
}

#' Test if in set
#' 
#' @export
#' @param y random vector, regular atomic vector, or \code{NULL}: the
#' values to be matched against.
#' @inheritParams rvmatch
"%*in*%" <- function (x, y) {
  if (! is.rv(x) && ! is.rv(y)) {
    return(.Primitive("%in%")(x, y))
  }
  z <- rvmatch(x, table=y, nomatch=0L)
  return(z > 0)
}
