# ========================================================================
# cbind  -  column bind for rvs
# ========================================================================
# Note: It's inconvenient that we cannot call the generic (default) cbind
#       if class attributes are set.
#

# DEBUG: cbind(1, rvnorm(1), rvnorm(1)) causes an error msg
#   In cbind(v, unclass(x[[i]]), deparse.level = deparse.level) :
#     number of rows of result is not a multiple of vector length (arg 2)



#' Combine random vectors by columns or rows
#' 
#' Combines random vectors by columns (\code{cbind.rv}) or rows
#' (\code{rbind.rv}).
#' 
#' See \link{cbind} and \link{rbind} for details.
#' 
#' @aliases cbind.rv rbind.rv
#' @param \dots vectors or matrices, can be rv objects
#' @param deparse.level (passed on to cbind)
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
#'   y <- rvnorm(10)
#'   cbind.rv(x, y)
#'   rbind.rv(x, y)
#' 
#' @method cbind rv
cbind.rv <- function(..., deparse.level = 1)
{
  if (deparse.level != 1) 
    .NotYetUsed("deparse.level != 1")
  x <- list(...)
  if (length(x)<1) return(NULL)
  v <- NULL
  for (i in seq(along=x)) {
    v <- cbind(v, unclass(x[[i]]), deparse.level=deparse.level)
  }
  class(v) <- class(rv())
  return(v)
}


## ========================================================================
## rvbind.rv  -  row bind for rvs
## ========================================================================
##

rbind.rv <- function(..., deparse.level = 1)
{
  if (deparse.level != 1) 
    .NotYetUsed("deparse.level != 1")
  x <- list(...)
  if (length(x)<1) return(NULL)
  v <- NULL
  for (i in seq(along=x)) {
    v <- rbind(v, unclass(x[[i]]), deparse.level=deparse.level)
  }
  class(v) <- class(rv())
  return(v)
}

