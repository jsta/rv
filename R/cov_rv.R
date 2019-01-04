# ========================================================================
# cov  -  short description
# ========================================================================

#' Calculate the covariance of an rv object
#' 
#' @importFrom stats cov
#' @rdname distrib_rv
#' @inheritParams stats::cov
#' @param \dots arguments passed to stats::cov
#' @export
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
