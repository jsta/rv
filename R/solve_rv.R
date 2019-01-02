#' Random Vectors
#' 
#' \code{solve.rv}
#' 
#' \code{solve.rv} is the rv-object compatible version of the function
#' \code{solve}.
#' 
#' For details of the function, see \code{\link{solve}}.
#' 
#' @param a a square random vector containing the coefficients of the linear
#' system
#' @param b a square random vector giving the right-hand side(s) of the linear
#' system
#' @param \dots further arguments passed to \code{solve}
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{solve}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @export
solve.rv <- function (a, b, ...) {
  rvmapply(base::solve, a=a, b=b, ...)
}
