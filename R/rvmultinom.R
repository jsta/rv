#' Generate Random Variables from a Multinomial Sampling Model
#' 
#' Generates a random vector from a multinomial sampling model.
#' 
#' The length of \code{prob} determines the number of bins.
#' 
#' The vector \code{prob} will be normalized to have sum 1.
#' 
#' If \code{length(prob)} is two, \code{rvbinom} is called instead.
#' 
#' NOTE. Case of random \code{n} or \code{size} or \code{prob} --- not yet
#' optimized for speed.
#' 
#' @param n integer, number of random variables to generate
#' @param size integer or integer-valued rv: the number of trials (size of each
#' sample)
#' @param prob vector (of length at least 3) prior probabilities of successes
#' of each trial (may be constant or an rv object)
#' @return A random array of dimensions \code{length(prob)} times \code{n}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   y <- rvmultinom(n=3, size=1, prob=c(0.20, 0.30, 0.50))
#' 
#' @export rvmultinom
rvmultinom <- function(n=1, size=1, prob) {
  if (length(prob)<=1) {
    return(rvbinom(n=n, size=size, prob=prob))
  }
  if (anyisrv(n, size, prob)) {
    r <- rvmapply(rmultinom, n=n, size=size, prob=prob)
  } else {
    n.sims <- getnsims()
    s <- rmultinom(n=n*n.sims, size=size, prob=prob)
    dim(s) <- c(length(s) %/% n.sims, n.sims)
    r <- rvsims(t(s))
    dim(r) <- c(length(prob), n)
    dimnames(r) <- list(names(prob), NULL)
  }
  return(r)
}

