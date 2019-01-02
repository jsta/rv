# rv-Math.R - standard math functions for the rv class



#' Mathematical functions and Operators for rv Objects
#' 
#' Mathematical functions and operators adapted to work with random variable
#' (rv) objects.
#' 
#' The operator method preserves the names of the longer vector (or those of
#' the first if the lengths match).
#' 
#' @aliases Math.rv Ops.rv !.rv Math.rvsim Ops.rvsim cumsum.rv cumprod.rv
#' cummin.rv cummax.rv
#' @param x object
#' @param \dots further arguments passed to or from other methods
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
#'   -x
#'   names(x) <- paste("x[", seq_along(x), "]", sep="")
#'   x + 1:10
#'   1:2 + x
#'   cumsum(x)
#'   cumprod(exp(x))
#' 
#' @method Math rv
Math.rv <- function(x, ...) {
  # Componentwise operation
  X <- x # Preserve class and other attributes
  for (i in seq_along(x)) {
    x <- X[[i]] 
    X[[i]] <- NextMethod()
  }
  return(X)
}

# cumsum, cumprod, cummax, cummin
#' @method cumsum rv
cumsum.rv <- function (x)
{
  simapply(x, cumsum)
}

#' @method cumprod rv
cumprod.rv <- function (x)
{
  simapply(x, cumprod)
}

#' @method cummin rv
cummin.rv <- function (x)
{
  simapply(x, cummin)
}

#' @method cummax rv
cummax.rv <- function (x)
{
  simapply(x, cummax)
}


# ----------------
# end of rv-Math.R
# ----------------
