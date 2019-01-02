#' Generate a Random Vector from an Empirical Distribution
#' 
#' \code{rvboot} generates a random vector of the same length as data from the
#' empirical distribution of the data.
#' 
#' \code{rvboot}
#' 
#' @param data A vector of constants
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   y <- rnorm(30) # Some data: 30 draws from standard normal.
#'   x <- rvboot(y) # A random vector of length 30 (each component has the same distribution)
#'   print(mean(x)) # Bootstrap estimate of the mean.
#'   print(sd.rv(x))   # Bootstrap estimate of the sd.
#'   rvinci(mean(x), 0) # Hypothesis test: mean of x is zero (at 5% level) FALSE => reject.
#' 
#' @export rvboot
rvboot <- function (data) {
#  empirical (bootstrap) distribution
#
  n.sims <- getnsims()
  n <- n.sims*length(data)
  s <- matrix(sample(data, size=n, replace=TRUE), nrow=n.sims)
  r <- rvsims(s)
  dim(r) <- dim(data)
  return(r)
}

