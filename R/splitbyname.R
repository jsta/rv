#' Split a vector based on the names of the components
#' 
#' \code{splitbyname} is a utility function that splits the given vector based
#' on the names of the components and returns a named list of arrays and
#' vectors.
#' 
#' The names are supposed to be of the format `name[index]`, for example
#' `alpha[1,1]`, `beta[1]`, etc.
#' 
#' A name without brackets is equivalent to a name with `[1]`.
#' 
#' The dimension attribute will not be set in case of vectors.
#' 
#' @param x a vector or a list with the name attributes set
#' @return A list of arrays and vectors.  Missing entries in the arrays and
#' vectors are filled in with \code{NA}s.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @keywords manip
#' @examples
#' 
#'   x <- structure(c(1,3), names=c("x[1,1]", "x[3,3]"))
#'   splitbyname(x) # yields a list containing a 3x3 matrix
#' 
#' @export splitbyname
splitbyname <- function (x) {
  a <- split(x, f = .shortnames(x))
  lapply(a, .setDimensionByName)
}

