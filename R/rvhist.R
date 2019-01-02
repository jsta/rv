# ========================================================================
# rvhist  -  plot histograms of the simulations of the random components
# ========================================================================
#



#' Histogram of Distributions of Components of a Random Vector
#' 
#' \code{rvhist} shows a grid of histograms of simulations of the components of
#' a random vector.
#' 
#' Outputs a histogram using the \code{hist} function with the option
#' \code{freq=FALSE}. This can be overridden by specifying the argument
#' \code{freq} or \code{prob}.  See the function \code{hist} for details.
#' 
#' @param x an rv object
#' @param \dots further arguments passed to the function \code{hist}
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @export rvhist
#' @importFrom graphics par
rvhist <- function (x, ...)  {
    if (!is.null(dim(x))) 
        par(mfcol = dim(x))
    mfcol <- par("mfcol")
    n <- prod(mfcol)
    n <- min(n, length(x))
    a <- list(...)
    if (is.null(a$freq) && is.null(a$prob)) {
      a$freq <- FALSE
    }
    if (is.null(a$breaks)) {
      a$breaks <- "fd"
    }
    make.main <- is.null(a$main)
    make.xlab <- is.null(a$xlab)
    lab <- deparse(substitute(x))
    x.names <- paste(lab, .dimindex(x), sep = "")
    for (i in 1:n) {
        a$x <- sims(x[i])
        if (make.main || make.xlab) {
            this.name <- x.names[i]
            if (make.xlab) {
                a$xlab <- this.name
            }
            if (make.main) {
                a$main <- paste("Histogram of", this.name)
            }
        }
        do.call("hist", a)
    }
}


