#' Sample from an arbitrary density function using grid approximation
#' 
#' \code{rvdens} generates a random vector where each simulation comes from a
#' Bernoulli sampling distribution.
#' 
#' 
#' @param n number of random scalars to draw
#' @param FUN density function
#' @param range range to discretize over
#' @param unitprecision how many points per unit length
#' @param \dots other arguments passed on to \code{FUN}
#' @return A random vector (an rv object) of length \code{n}.
#' @note The resulting vector will not be independent and identically
#' distributed Bernoulli unless \code{prob} is a fixed number.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvdens(FUN=stats:::dnorm, range=c(-5, 5), unitprecision=10)
#'   y <- rvnorm(1) ## Should be close to x
#' 
#' @export rvdens
rvdens <- function(n=1, FUN, range, unitprecision=10, ...) {
  # NAME
  #   rvdensity - Sample from a given univariate density using a grid approximation
  # ARGUMENTS
  #   n : number of independent random vector components to draw
  #   FUN : density function, must be vectorized
  #   range : range for the grid
  #   unitprecision : number of points per unit
  #   ... : other arguments passed to [FUN].
  #   
  grid <- seq(from=range[1], to=range[2], by=1/unitprecision)
  prob <- FUN(grid, ...)
  n.sims <- getnsims()
  s <- sample(grid, size=n*n.sims, prob=prob, replace=TRUE)
  noise <- runif(n*n.sims, -0.5/unitprecision, 0.5/unitprecision)
  rvsims(matrix(s+noise, nrow=n.sims))
}


