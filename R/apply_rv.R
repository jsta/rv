#' Apply Functions over Margins of Random Arrays
#' 
#' The \code{rv}-compatible version of \code{apply}
#' 
#' This is the rv-compatible version of the function \code{\link{apply}}.
#' 
#' @param X a random array
#' @param MARGIN subscripts.
#' @param FUN function.
#' @param \dots optional arguments to \code{FUN}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{apply}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords manip
#' @examples
#' 
#' \dontrun{
#'   x <- rvmatrix(rvnorm(12), nrow=3, ncol=4)
#'   print(apply.rv(x, 1, sum))
#' }
#' 
#' @export apply.rv
apply.rv <- function (X, MARGIN, FUN, ...) { ## CHECK
  ## NAME
  ##   apply.rv - Apply Functions Over Random Array Margins
  if (length(dim(X)) == 0) {
    stop("dim(X) must have a positive length")
  }
  a <- apply(X, MARGIN, FUN, ...)
  cat("NOTE: 'apply.rv' does not yet set the dimensions of the result (TODO)")
  return(a)
}

