# ========================================================================
# hist.rv  -  histogram, adapted for rv's
# ========================================================================



#' Histogram of a random vector
#' 
#' \code{hist.rv} shows a grid of histograms generated from random draws of the
#' random vector argument.
#' 
#' 
#' @param x an object
#' @param grid a vector of two numbers, indicating the size of the grid to plot
#' the histograms
#' @param xlim x limits
#' @param main main title
#' @param freq logical; if \code{FALSE}, plots as probability density, as it
#' should.
#' @param \dots Other arguments passed on to \link{hist}
#' 
#' @importFrom graphics hist
#' 
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   \dontshow{require(rv)}
#'   x <- rvnorm(30)
#'   hist(x)
#' 
hist.rv <- function(x, grid=c(4,5), xlim=x.range, main=paste(xname,"simulation"), freq=FALSE, ...) {
  par(mfrow=grid)
  #  l <- length(x)
  #  par(mfrow=c(l %% 3 + l %/% 3, 3))
  xname <- deparse(substitute(x))
  grid.n <- grid[1]*grid[2]
  s <- sims(x)
  x.range <- c(min(s),max(s))
  if (grid.n<1)
    stop("Bad grid")
  s <- s[1:grid.n,]
  for (i in 1:grid.n) {
    # truehist(s[i,], xlim=xlim, main, ...)
    hist(s[i,], xlim=xlim, main=main, freq=freq, ...)
  }
}

