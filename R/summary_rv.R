#' Numerical Summary of a Random Variable
#' 
#' Gives a numerical summary of the random variable in the format of a data
#' frame.
#' 
#' The objects are first coerced to \code{rvsummary} objects, then passed on to
#' the \code{summary.rvsummary} method, which creates a nicely formatted data
#' frame of the object.
#' 
#' @aliases summary.rv summary.rvfactor summary.rvmixed summary.rvsummary
#' summary.rvsummary_numeric summary.rvsummary_integer
#' summary.rvsummary_logical summary.rvsummary_rvfactor
#' @param object object to summarize
#' @param all.levels show summary for all levels even if there are too many to
#' display in one line
#' @param \dots rv object
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords internal
#' @examples
#' 
#'   x <- rvarray(rvnorm(6), c(2,3))
#'   summary(x)
#'   summary(as.rvsummary(x))
#'   summary(rvfactor(trunc(x)))
#' 
#' @method summary rv
summary.rv <- function (object, ...) {
  summary(as.rvsummary(object), ...)
}

#' @method summary rvfactor
#' @export
summary.rvfactor <- function (object, all.levels=TRUE, ...) {
  summary(as.rvsummary(object), all.levels=all.levels, ...)
}
