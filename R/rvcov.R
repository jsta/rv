# ========================================================================
# rvcov  -  covariance matrix
# ========================================================================
#



#' Covariance Between Components of Random Vectors
#' 
#' \code{rvcov}
#' 
#' \code{rvcov}
#' 
#' @param x a random vector
#' @param y (optional) a random vector
#' @param \dots further arguments passed to or from other methods
#' @return A covariance matrix.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(mean=1:3)
#'   y <- rvnorm(mean=2:4)
#'   rvcov(x,y)
#'   rvcov(x,x)
#' 
#' @export rvcov
rvcov <- function(x, y=NULL, ...)
{
  if (is.null(y)) {
    cov(sims(x), y, ...)
  } else {
    cov(sims(x), sims(y), ...)
  }
}

