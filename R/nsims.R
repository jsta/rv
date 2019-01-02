#' Number of simulations stored in each component of an rv object
#' 
#' \code{rvnsims} returns the number of simulations stored in each component of
#' its argument; \code{setnsims} sets the default number of simulations;
#' \code{getnsims} retrieves the default number of simulations.
#' 
#' If the argument is a non-rv numeric vector, \code{rvnsims} returns 1
#' (corresponding to a `constant') for each component.
#' 
#' The minimum number of default simulations is 2.
#' 
#' @aliases rvnsims rvnsims.rv rvnsims.rvsummary setnsims getnsims
#' @param x an rv object.
#' @param n.sims default number of simulations; must be at least 2.
#' @return \code{rvnsims}: a vector of integers.
#' 
#' \code{setnsims}: \emph{previously set} default number of simulations.
#' 
#' \code{getnsims}: (integer) currently set default number of simulations.
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
#'   rvnsims(1.23)             # 1
#'   x <- rvnorm(1)            # an rv
#'   rvnsims(x)                # equal to setnsims()
#'   rvnsims(x)==nrow(sims(x)) # TRUE
#'   rvnsims(x)==getnsims()    # TRUE
#'   setnsims(1000)            # set n.sims to 1000
#'   n.sims <- setnsims(10000) # s is now 1000
#'   print(getnsims())         # prints 10000
#'   setnsims(n.sims)          # restore the number of simulations back to 1000
#' 
#' @export rvnsims
rvnsims <- function (x) {
  UseMethod("rvnsims")
}

rvnsims.rv <- function (x) {
  sapply(unclass(x), length)
}

rvnsims.rvsummary <- function (x) {
  unlist(rvattr(x, "n.sims"), use.names=TRUE)
}

rvnsims.default <- function (x) {
  if (!(is.atomic(x) || is.recursive(x))) {
    stop("rvnsims: no applicable method for class '", class(x), "'")
  }
  rep.int(1, length(x))
}

#' @export
setnsims <- function (n.sims) {
  ## setnsims - get or set the default number of simulations (a global variable)
  if (length(n.sims)>0 && is.numeric(n.sims) && (!is.na(n.sims[1])) && n.sims[1]>=2) {
    n.sims <- as.integer(ceiling(n.sims[1]))
    oldn.sims <- rvpar("n.sims")
    rvpar(n.sims=n.sims)
  } else {
    stop('Invalid number of simulations (must be at least 2)', n.sims[1])
  }
  return(oldn.sims)
}

getnsims <- function () {
  n.sims <- rvpar("n.sims")
  if (!is.integer(n.sims) || n.sims<2) {
    stop("Invalid number of simulations - please set with setnsims(...)")
  }
  return(n.sims)
}

