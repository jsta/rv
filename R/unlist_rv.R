#' Flatten Lists Containing rv Objects
#' 
#' Given a list structure \code{x}, \code{unlist} simplifies it to produce a
#' vector which contains all the atomic components (\emph{containing rv
#' objects}) which occur in \code{x}.
#' 
#' This is the rv-compatible version of the function \code{\link{unlist}}.
#' 
#' Since \code{unlist} is not a generic function, the whole name
#' \code{unlistrv} must be specified when calling the function when \code{x} is
#' an 'rv' object.
#' 
#' @param x An R object, typically a list or vector (containing rv objects)
#' @param recursive logical. Should unlisting be applied to list components of
#' x?
#' @param use.names logical. Should names be preserved? (now fixed to TRUE)
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{unlist}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords manip
#' @examples
#' 
#'   x <- list(a=rvnorm(2), b=rvnorm(3))
#'   print(unlistrv(x))
#' 
#' @export unlistrv
unlistrv <- function (x, recursive = TRUE, use.names = TRUE) {
  # Flatten Lists That (May) Contain Random Variables 
  y <- NULL
  ix <- seq(along=x)
  xn <- names(x)
  .paste <- function (name, x) {
     nbrs <- .dim.index(x, leftadjust=FALSE)
     paste(name, nbrs, names(x), sep="")
  }
  if (recursive) {
    for (i in ix) {
      nx <- xn[i]
      if (use.names && is.null(nx)) nx <- "."
      if (!is.rv(x[[i]]) && is.list(x[[i]])) {
        new.y <- unlistrv(x[[i]], recursive=TRUE, use.names=use.names)
      } else {
        new.y <- x[[i]]
      }
      if (is.null(names(new.y))) {
        new.names <- .paste(nx, new.y)
      } else {
        new.names <- paste(nx, ".", names(new.y), sep="")
      }
      yn <- names(y)
      y <- c(y, new.y)
      names(y) <- c(yn, new.names)
    }
  } else {
    for (i in ix) {
      nx <- xn[i]
      new.y <- x[[i]]
      new.names <- .paste(nx, new.y)
      yn <- names(y)
      y <- c(y, new.y)
      names(y) <- c(yn, new.names)
    }
  }
  return(y)
}
