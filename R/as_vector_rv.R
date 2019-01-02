#' Coerce an rv object
#' 
#' \code{as.vector.rv} coerces a given \code{rv} object into a vector; matrices
#' lose their dimension attributes, but \code{rv} objects stay as \code{rv}
#' objects (since they are considered to be ``vectors'').
#' 
#' \code{as.vector.rv} removes the dimension attribute and returns the rv
#' object.  Needed for compatibility with code that uses \code{as.vector}.
#' 
#' @param x an object
#' @param mode (currently not used)
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvmatrix(rvnorm(10), 2, 5)
#'   as.vector(x)
#' 
as.vector.rv <- function(x, mode="any") {
  a <- attributes(x) 
  x <- lapply(unclass(x), as.vector, mode=mode)
  attributes(x) <- a
  dim(x) <- NULL
  return(x)
}


