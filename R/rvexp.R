#' Generate Random Vectors from an Exponential Sampling Model
#' 
#' \code{rvexp}
#' 
#' \code{rvexp}
#' 
#' @param n integer: number of variables to generate
#' @param rate prior distribution for the rate parameter (constant or random)
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   y <- rvexp(1, rate=rvexp(1)) # What marginal distribution does y have now?
#' 
#' @export rvexp
rvexp <- function (n=1, rate=1) {
  rvvapply(rexp, n.=n, rate=rate)
}


