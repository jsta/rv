#' Generate Random Variables from a Binomial Sampling Model
#' 
#' Generates a random vector from a binomial sampling model.
#' 
#' \code{rvbinom} generates a random vector with given length, the distribution
#' for size and the distribution for the probability of success.
#' 
#' @param n integer, number of random variables to generate
#' @param size integer or integer-valued rv: the number of trials (size of each
#' sample)
#' @param prob prior probability of success of each trial (may be constant or
#' an rv object)
#' @return An rv object.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   s <- 1+rvpois(1,lambda=3)        # A prior distribution on the 'size' parameter.
#'   rvbinom(1, size=s, prob=0.5)     # The 'size' is random.
#'   p <- rvbinom(1, 10, prob=0.5)/10 # Prior probability of success.
#'   rvbinom(1, size=10, prob=p)      # Now the probability is random.
#'   rvbinom(1, size=s, prob=p)       # Both the size and the probability are random.
#' 
#' @export rvbinom
rvbinom <- function (n=1, size, prob) {
  rvvapply(rbinom, n.=n, size=size, prob=prob)
}

