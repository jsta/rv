#' Replicate Elements of Random Vectors
#' 
#' Transpose a random array by permuting its dimensions and optionally resizing
#' it.
#' 
#' This is the rv-compatible version of the function \code{\link{rep}}.
#' 
#' Since \code{rep} is not a generic function, the whole name \code{rep.rv}
#' must be specified when calling the function when \code{x} is an 'rv' object.
#' 
#' @param x a random vector to be replicated
#' @param times number of replications
#' @param \dots further arguments passed to \code{rep}
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{rep}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords manip
#' @examples
#' 
#'   print(rep(rvnorm(1), times=4))
#' 
rep.rv <- function (x, times, ...) {
  if (! is.rv(x)) return(rep(x, times, ...))
  if (missing(times)) {
    a <- rep(unclass(x), ...)
  } else {
    a <- rep(unclass(x), times, ...)
  }
  class(a) <- class(x)
  return(a)
}
