#' Categorical Random Variables (Random Factors)
#' 
#' Creates or tests for objects of type ``\code{rvfactor}''.
#' 
#' Internally random factors are integer-valued just like regular factors in R.
#' 
#' The number of levels to print when \code{all.levels==FALSE} can be set by
#' \code{rvpar(max.levels=...)}. By default this is set to 10.
#' 
#' @aliases rvfactor rvfactor.rv is.rvfactor as.rvfactor as.rv.rvfactor
#' print.rvfactor
#' @param x object to be coerced or tested.
#' @param all.levels logical; whether to print all levels or not (see below for
#' details)
#' @param \dots other arguments
#' @return \code{rvfactor}: an \code{rvfactor} object.
#' 
#' \code{is.rvfactor}: \code{TRUE} or \code{FALSE}.
#' 
#' \code{as.rv.rvfactor}: an \code{rv} object.
#' 
#' \code{as.rvfactor.rv}: an \code{rvfactor} object.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   # Probabilities of each integer of trunc(Z) where Z ~ N(0,1) ?
#'   x <- rvnorm(1)
#'   rvfactor(trunc(x))
#'   rvfactor(x>0)
#'   rvfactor(rvpois(1, lambda=0.5))
#' 
#' @export rvfactor
rvfactor <- function (x, ...) {
  UseMethod("rvfactor")
}

#' @rdname rvfactor
#' @param levels factor levels (labels for the levels)
#' @method rvfactor default
rvfactor.default <- function (x, levels=NULL, ...) {
  f <- as.factor(x)
  a <- sims(as.rv(as.integer(f)))
  rvf <- rvsims(a)
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- if (is.null(levels)) {
    levels(f)
  } else {
    levels
  }
  return(rvf)
}

#' @method rvfactor rv
rvfactor.rv <- function (x, levels=NULL, ...) {
  a <- sims(x)
  f <- as.factor(a)
  levs <- levels(f)
  rvf <- rvsims(array(as.integer(f), dim(a)))
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- if (is.null(levels)) {
    levs
  } else {
    levels
  }
  return(rvf)
}

is.numeric.rvfactor <- function (x) {
  FALSE
}

is.rvfactor <- function (x) {
  UseMethod("is.rvfactor")
}

#' @method is.rvfactor rvfactor
is.rvfactor.rvfactor <- function (x) {
  TRUE
} 

is.rvfactor.rv <- function (x) {
  all(rvsimapply(x, is.factor))
} 

is.rvfactor.default <- function (x) {
  FALSE
} 

as.rvfactor <- function (x, ...)
{
  if (is.rvfactor(x)) x else rvfactor(x)
} 


as.rv.rvfactor <- function (x, ...) {
  return(x)
  attr(x, "levels") <- NULL
  clx <- class(x)
  clx <- clx[clx!="rvfactor"]
  class(x) <- clx
  return(x)
}

"[.rvfactor" <- function (x, ..., drop = FALSE) {
  y <- NextMethod("[")
  attr(y, "levels") <- attr(x, "levels")
  class(y) <- oldClass(x)
  lev <- levels(x)
  if (drop) {
    exclude <- if (any(is.na(levels(x)))) { NULL } else { NA }
    y <- factor(y, exclude=exclude)
  }
  return(y)
}
