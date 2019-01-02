#' Distribution of the Arithmetic Mean of a Random Vector
#' 
#' \code{mean.rv} computes the distribution of the arithmetic average of its
#' argument \code{rv} object.
#' 
#' \code{mean} gives the distribution (that is, a random variable object) of
#' the statistic \eqn{\frac{1}{n}\sum_{i=1}^n x_i}{sum x_i/n}
#' (\code{sum(x)/length(x)}).
#' 
#' In particular, \code{mean(x)} of a random vector \code{x} of length one is
#' equal to \code{x} as it would be in the case of numerical x.
#' 
#' To find the expectation of a random vector \code{x} (that is, the individual
#' means of random components in a vector), use \code{rvmean(x)} (same as
#' \code{E(x)} and \code{Pr(x)}).
#' 
#' @param x an object
#' @param \dots further arguments passed to or from other methods
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{rvmean}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   y <- rvnorm(10, mean=0, sd=1)
#'   m1 <- mean(y)
#'   m2 <- rvnorm(1, mean=0, sd=1/sqrt(10))
#'   print(c(m1, m2)) # should have the same distribution
#' 
mean.rv <- function(x, ...) {
  rvsims(rowMeans(sims(x), ...))
}


