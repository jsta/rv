#' Generate Random Vectors from a Beta Sampling Model
#' 
#' \code{rvbeta} generates a random vector from the beta sampling model;
#' 
#' \code{rvnbeta(n, a, b)} ("neutral" Beta distribution) is equivalent to
#' \code{rvbeta(n, 1/3+a, 1/3+b)}.
#' 
#' 
#' @aliases rvbeta rvnbeta
#' @param n integer, number of random variables to generate
#' @param shape1 positive number or rv, 1st shape parameter
#' @param shape2 positive number or rv, 2nd shape parameter
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'    n <- 12         # sample size
#'    y <- (0:(n-1))  # observations
#'    a <- b <- 1/3   # the neutral beta prior
#'    rvbeta(1, shape1=a+y, shape2=b+n-y)
#'    rvnbeta(1, shape1=y, shape2=n-y)
#' 
#' @export rvbeta
rvbeta <- function (n=1, shape1, shape2) {
  rvvapply(rbeta, n.=n, shape1=shape1, shape2=shape2)
}

