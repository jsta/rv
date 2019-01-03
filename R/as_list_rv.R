#' Coerce a random vector object to a list
#' 
#' \code{as.list.rv} coerces a given \code{rv} object into a list.
#' 
#' Each component of the argument is extracted into a component of an enclosing
#' list, which is returned.
#' 
#' @param x an rv object
#' @param \dots arguments passed on to other methods
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(10)
#'   L <- as.list(x)
#' 
#' @export
#' @method as.list rv
as.list.rv <- function (x, ...) { ## 
  L <- vector(mode="list", length=length(x))
  for (i in seq_along(x)) {
    L[[i]] <- x[i]
  }
  return(L)
}
