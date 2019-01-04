#' Random Matrix Multiplication
#' 
#' Multiplies two random matrices, if they are conformable.  If one argument is
#' a vector, it will be coerced to either a row or column matrix to make the
#' two arguments conformable.  If both are vectors it will return the inner
#' product.
#' 
#' Optimized internally for the case of random matrix multiplied by a constant
#' column vector.
#' 
#' @param x,y numeric or complex matrices or vectors.
#' @name matmult
#' @return The (distribution of the) matrix product.  Use \code{\link{drop}} to
#' get rid of dimensions which have only one level.
#' @seealso \code{\link{matrix}}, \code{\link{Ops}}, \code{\link{diag}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords array arith
#' @examples
#' 
#' x <- 1:4
#' (z <- x %*% x)    # scalar ("inner") product (1 x 1 matrix)
#' drop(z)             # as scalar
#' 
#' y <- diag(x)
#' z <- matrix(1:12, ncol = 3, nrow = 4)
#' y %*% z
#' y %*% x
#' x %*% z
#' @export
'%*%.rv' <- function(x, y) { ## CHECK: TTRY if is.constant(x) => normal
  # ========================================================================
  # %*% - matrix product
  # ========================================================================
  # DEBUG: how to make this work with outer()?
  return("%**%"(x, y))
}

#' Constant matrix times a rv vector
#' 
#' @inheritParams base::`%*%`
#' @rdname matmult
#' @export
`%**%` <- function(x, y) { ## CHECK: TTRY if is.constant(x) => normal
  if (! is.rv(x) && ! is.rv(y)) {
    return(.Primitive("%*%")(x, y))
  }
  d <- dim(y)
  if (! is.rv(x) && (is.null(d)) || (length(d) == 2 && d[2] == 1)) {
    n.sims <- .Primitive("max")(rvnsims(x), rvnsims(y), na.rm=FALSE)
    ysim <- sims(as.rv(y), dimensions=TRUE, n.sims=n.sims)
    # Typical case: constant matrix times a rv vector
    AB <- t(.Primitive("%*%")(x, t(ysim)))
    rvsims(AB)
  } else {
    rvmapply(base::crossprod, x=t(as.rv(x)), y=as.rv(y))
  }
}
