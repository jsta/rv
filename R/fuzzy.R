#' Fuzziness
#' 
#' Tests whether an object is "fuzzy", i.e.  a logical random scalar that has
#' probability strictly between zero and one (not strictly true nor strictly
#' false).
#' 
#' 
#' @aliases fuzzy is.fuzzy is.fuzzy.rv
#' @param x an object, random or constant
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- as.logical(rvbern(1,0.4)) # a logical random variable
#'   is.fuzzy(x) # TRUE, since x is logical and not constant
#'   is.fuzzy(x<2) # FALSE, since x is less than 2 with probability one
#'   is.fuzzy(rvnorm(1)) # FALSE, since it's not a probability
#'   is.fuzzy(TRUE) # FALSE, since TRUE is strictly TRUE
#'   is.fuzzy(1) # FALSE, since 1 is not a logical variable
#' 
#' @export
is.fuzzy <- function (x) {
  UseMethod("is.fuzzy")
}

#' @method is.fuzzy rv
#' @export
is.fuzzy.rv <- function (x) {
  # NAME
  #  is.fuzzy - Is a Vector Component Logical But Random
  # 
  component.is.logical <- rvsimapply(x, is.logical)
  component.prop <- rvmean(x)
  (component.is.logical & component.prop>0 & component.prop<1)
}

#' @method is.fuzzy default
#' @export
is.fuzzy.default <- function (x)
{
  return(FALSE)
}
