# ========================================================================
# cov  -  short description
# ========================================================================

#' @importFrom stats cov
cov.rv <- function(x, y=NULL, ...)  ## EXPORT cov.rv
{
  if (!is.matrix(x)) {
    if (is.null(y)) {
      stop("supply both x and y or a matrix-like x")
    }
    x <- as.vector(x)
  }
  if (is.rvobj(x) || is.rvobj(y)) {
    rvmapply(cov, x=x, y=y, ...)
  } else {
    cov(x=x, y=y, ...)
  }
}



