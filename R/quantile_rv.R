# ========================================================================
# quantiles of a random vector
# ========================================================================
#



#' Distribution of a Quantile of a Random Vector
#' 
#' \code{quantile.rv} returns the distribution of the quantile of a random
#' vector (as a random variable).
#' 
#' 
#' @param x an object
#' @param \dots further arguments passed to or from other methods
#' @return A random vector (rv object) with components giving the distribution
#' of the desired quantiles.
#' @note \code{quantile.rv} does not return the simulated quantiles of the
#' quantiles of the argument \code{x}.  This is done by
#' \code{\link{rvquantile}}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(30)
#'   quantile(x)
#' 
#' @method quantile rv
#' @export
quantile.rv <- function(x, ...) {
  simapply(x, quantile, ...)
}

