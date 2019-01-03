#' Generate Random Vectors from a Discrete Sampling Model
#' 
#' Generates random variables from a discrete distribution (from a finite
#' population with replacement).
#' 
#' Computes a random vector of length \code{n}, consisting of identicallly
#' distributed discrete random scalars with the discrete distribution with
#' values \code{x} and corresponding probabilities \code{prob}.  If \code{prob}
#' is not given, all values are considered equally distributed.
#' 
#' @param n integer: number of scalars to generate
#' @param x values of the distribution
#' @param prob probabilities (optional, default: all equal)
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples \dontrun{
#' 
#'   # 8 people draw a number each from 1..10 with replacement.
#'   # What is the probability that the highest number of the eight is "10"?
#'   u <- rvdiscrete(n=8, x=1:10) # 8 iid variables from the discrete uniform 1:10.
#'   Pr(max(u)==10)
#'   # What is the probability that the person with the 3rd smallest number
#'   # has at least "3"?
#'   s <- sort(u) # order distribution
#'   Pr(s[3]>=3)
#'   }
#' 
#' @export rvdiscrete
rvdiscrete <- function (n=1, x, prob=NULL) {
  n.sims <- getnsims()
  rvsims(matrix(sample(x=x, size=n * n.sims, prob=prob, replace=TRUE), nrow=n.sims))
}

