#' Convert Numeric to Random Factor
#' 
#' Convert implements the `cut' function using random variables.
#' 
#' 
#' @aliases rvcut rvcut.rv
#' @param x a plain or a random vector which is to be converted to a factor by
#' cutting.
#' @param \dots arguments passed to the function \code{\link{cut}}.
#' @return A random factor.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{rvfactor}}, \code{\link{cut}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   rvcut(rvnorm(1), breaks=c(-Inf,-2,-1,0,1,2,Inf))
#' 
#' @export rvcut
rvcut <- function (x, ...) {
  UseMethod("rvcut")
}

#' @method rvcut default
#' @export
rvcut.default <- function (x, ...) {
  f <- cut(x, ...)
  levs <- levels(f)
  rvf <- rvsims(as.integer(f))
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- levs
  return(rvf)
}

#' @method rvcut rv
#' @export
rvcut.rv <- function (x, ...) {
  a <- sims(x)
  f <- cut(a, ...)
  levs <- levels(f)
  rvf <- rvsims(array(as.integer(f), dim(a)))
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- levs
  return(rvf)
}


