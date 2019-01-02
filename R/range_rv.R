#' Distribution of the Range of a Random Vector
#' 
#' \code{range.rv} returns a 2-component random vector containing the
#' distributions of the minimum and the maximum values of all the given
#' arguments.
#' 
#' This is the rv-compatible version of the function \code{\link{range}}.
#' 
#' @param \dots further arguments passed to or from other methods
#' @param na.rm logical, indicating if \link{NA}s should be omitted
#' @param finite logical, indicating if all non-finite elements should be
#' omitted
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{quantile.rv}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(mean=1:10, sd=1)
#'   print(range(x))
#'   print(quantile(x, c(0,1)))
#' 
range.rv <- function(..., na.rm=FALSE, finite=FALSE) {
  sims <- sims(c(...)) # an L x n matrix of numbers
  m <- apply(sims, 1, 'range', na.rm=na.rm, finite=finite) # gives a 2xL matrix!!!
  r <- rvsims(t(m)) # Transpose to get an L x 2 matrix.
  names(r) <- c('min', 'max')
  return(r)
}
