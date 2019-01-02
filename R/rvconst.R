#' Random Vector with a Point-Mass Distribution
#' 
#' Coerces a given vector of constants into a random vector with 1 simulation
#' in each component.
#' 
#' Coerces a given vector of constants into a random vector with 1 simulation
#' in each component.
#' 
#' @param n integer: number of variables to generate
#' @param x a vector of constants
#' @return A random vector (rv object) of length \code{n}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvconst(x=1:3)
#'   c(x, 4)
#' 
#' @export rvconst
rvconst <- function(n=1, x=0) {
  lx <- length(x)
  v <- rv(lx)
  if (lx > 0) {
    v[1:lx] <- x # recycle
  }
  return(v)
}

