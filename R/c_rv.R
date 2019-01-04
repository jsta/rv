#' Combine values in an rv object
#' 
#' @description Concatentates random vectors.
#' 
#' @details NOTE: \code{recursive} has not yet been tested.
#' 
#' \code{cc} is a function that works for both non-rv and other vectors. To
#' make code compatible for both constant vectors and rv objects, one can use
#' \code{cc} instead of \code{c}.
#' 
#' @aliases c.rv c.rvsummary cc
#' @param \dots objects to be concatenated. Can be a mixture of constants and
#' rv objects.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(2)
#'   y <- rvbern(2, prob=0.5)
#'   z <- c(x, y)
#'   print(z)
#'   z1 <- cc(1, z)
#'   z2 <- c(as.rv(1), z)
#'   z3 <- c(as.rv(1), z)
#'   print(z1)
#'   print(z2)
#'   print(z3)
#' 
#' @inheritParams base::c
#' @export
cc <- function(..., recursive=FALSE) {
  args <- lapply(as.list(match.call())[-1], eval, envir=parent.frame())
  cls <- sapply(args, class)
  if ("rvsummary" %in% cls) {
    x <- c.rvsummary(..., recursive=recursive)
  } else if ("rv" %in% cls) {
    x <- c.rv(...)
  } else {
    x <- c(..., recursive=recursive)
  }
  return(x)
}

#' @inheritParams cc
#' @param recursive logical. If recursive = TRUE, the function recursively
#' descends through lists (and pairlists) combining all their elements into a
#' vector.
#' @rdname cc
#' @export
#' @method c rv
c.rv <- function(..., recursive=FALSE)
{
  ## a kludge to disable dispatching rv
  x <- c(list(NA), ..., recursive=recursive)[-1]
  class(x) <- class(rv())
  return(x)
}

c.rvsummary <- function(..., recursive=FALSE)
{
  ## a kludge to disable dispatching rv
  x <- c(list(NA), ..., recursive=recursive)[-1]
  class(x) <- "rvsummary"
  return(x)
}
