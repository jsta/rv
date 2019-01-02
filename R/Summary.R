#' Summary
#' 
#' Methods applied to rv objects for the generic functions \code{all},
#' \code{any}, \code{sum}, \code{min}, \code{max}.
#' 
#' 
#' @param \dots rv object
#' @param na.rm logical, remove NAs?
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords internal
#' @examples
#' 
#'   x <- rvnorm(5)
#'   all(x>0)
#'   any(x>2)
#'   sum(x)
#'   min(x)
#'   max(x)
#' @method Summary rv
Summary.rv <- function(..., na.rm=FALSE)
{
  S <- sims(c(...)) # an L x n matrix of numbers
  M <- apply(S, 1, .Generic, na.rm=na.rm)
  rvsims(M)
}
