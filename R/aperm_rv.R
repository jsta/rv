#' Random Array Transposition
#' 
#' Transpose a random array by permuting its dimensions and optionally resizing
#' it.
#' 
#' This is the rv-compatible version of the function \code{\link{aperm}}.  It
#' first applies
#' 
#' @param a the random matrix to be transposed
#' @param perm the subscript permutation vector. See the manual page for the
#' gneric method aperm.
#' @param \dots further arguments passed to \code{aperm}
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{aperm}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords manip
#' @examples
#' 
#'   x <- rvarray(rvnorm(24), dim=c(2,3,4))
#'   print(aperm(x))
#' 
aperm.rv <- function (a, perm, ...) {
  # NAME
  #   aperm.rv - Transpose a Random Array
  #
  A <- NextMethod()
  if (! is.rv(a)) {
    return(A)
  }
  class(A) <- class(a)
  return(A)
}


