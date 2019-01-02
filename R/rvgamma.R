#' Generate Random Variables from a Gamma Sampling Model
#' 
#' Generates random variables from a Gamma sampling model.
#' 
#' \code{rvngamma(n, shape, rate)} is equivalent to \code{rvgamma(n, 1/3 +
#' shape, rate)}.
#' 
#' @aliases rvgamma rvngamma
#' @param n integer: number of variables to generate
#' @param shape shape parameter, may be a rv
#' @param rate rate parameter, may be a rv
#' @param scale inverse of rate, may be specified optionally instead of rate
#' @return A random vector (rv object).
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   round(rvmedian(rvngamma(n=1, shape=1:10, rate=1)), 1) ## close to 1:10
#' 
#' @export rvgamma
#' @importFrom stats rgamma
rvgamma <- function (n=1, shape, rate = 1, scale = 1/rate)  {
  rvvapply(rgamma, n.=n, shape=shape, scale=scale)
}


