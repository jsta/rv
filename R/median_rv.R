# ========================================================================
# median.rv - median of a random vector
# ========================================================================



#' Distribution of the Sample Median
#' 
#' Compute the distribution sample median of the vector of values given as its
#' argument.
#' 
#' 
#' @param x a randomv vector containing the components whose distribution of
#' the median value is to be computed.
#' @param na.rm a logical value indicating whether \code{NA} values should be
#' stripped before the computation proceeds.
#' @param \dots further arguments passed to \code{median}
#' @seealso \code{\link{rvmedian}} for the componentwise medians.
#' \code{\link{quantile}} for general quantiles.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(10) ## A random vector of length 10.
#'   median(x)       ## A random scalar (vector of length 1).
#'   rvmedian(x)     ## A numeric vector of length 10.
#' 
#' @export
#' @importFrom stats median
#' @method median rv
median.rv <- function(x, na.rm=FALSE, ...) {
  simapply(x, stats::median, na.rm=na.rm, ...)
}
