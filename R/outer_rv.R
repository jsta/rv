#' Outer Product of Random Arrays
#' 
#' \code{outer.rv}
#' 
#' Implements the outer product for random arrays.
#' 
#' Note. \code{outer} is not a generic function; thus \code{outer(x)} will not
#' work if \code{x} is an rv object.  You must write \code{outer.rv(x)}
#' explicitly.
#' 
#' See the function \code{outer} for further details.
#' 
#' @param X First argument for function \code{FUN}
#' @param Y Second argument for function \code{FUN}; if missing, \code{X} is
#' used instead
#' @param FUN a function to use on the outer products; a character string or a
#' function
#' @param \dots optional arguments to be passed to \code{FUN}
#' @return A random array.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   #
#' 
#' @export outer.rv
outer.rv <- function (X, Y=NULL, FUN="*", ...) {
  # NAME
  #   outer.rv - 
  # 
  if (is.null(Y)) {
    rvmapply("outer", X, X, MoreArgs=list(FUN=FUN, ...))
  } else {
    rvmapply("outer", X, Y, MoreArgs=list(FUN=FUN, ...))
  }
}
