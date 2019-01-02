# ========================================================================
# abline.rv  -  a + bx
# ========================================================================
# change this to
#   mapply.rv("abline", ..., MoreArgs=NULL)
#



#' Add (Random) Straight Lines to a Plot
#' 
#' \code{abline.rv}, with random arguments (i.e. arguments of which at least
#' one is an \code{rv} object), plots a sample of lines corresponding to of
#' simulations of rv object \code{x}.  If the arguments are all numeric (none
#' is an \code{rv} object), the function call is passed on to \code{abline}.
#' 
#' This is a version of \code{abline} that accepts random variable objects for
#' the arguments \code{a}, \code{b}, \code{h}, or \code{v}.
#' 
#' The number of lines is determined by \code{rvpar("line.sample")}, default
#' 20.
#' 
#' See the original help page in package `graphics.'
#' 
#' @param a intercept
#' @param b slope
#' @param h y-value(s) horizontal line(s)
#' @param v x-value(s) horizontal line(s)
#' @param \dots further arguments passed to \code{\link{abline}}
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords aplot
#' @examples
#' 
#'   \dontrun{
#'      demo("rvexample1")
#'   }
#' 
#' @export abline.rv
abline.rv <- function (a = NULL, b = NULL, h = NULL, v = NULL, ...) { ## NEW
  if (! anyisrv(a, b, h, v)) {
    return(abline(a=a, b=b, h=h, v=v, ...))
  }
  line.sample <- rvpar("line.sample")
  if (!is.numeric(line.sample)) {
    stop("rvpar('line.sample') is not a number")
  }
  args <- list(FUN=abline, a=a, b=b, h=h, v=v)
  nulls <- sapply(args, is.null)
  nullArgs <- names(nulls)[nulls]
  MoreArgs <- list(...)
  args[nullArgs] <- NULL
  MoreArgs[nullArgs] <- list(NULL)
  args$SAMPLESIZE <- line.sample
  args$MoreArgs <- MoreArgs
  do.call(rvmapply, args=args)
  invisible()
}

