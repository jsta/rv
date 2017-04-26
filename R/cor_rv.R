# ========================================================================
# cor  -  correlation function
# ========================================================================

cor.rv <- function(x, y=NULL, ...)  {
  if (!is.matrix(x)) {
    if (is.null(y)) {
      stop("supply both x and y or a matrix-like x")
    }
    x <- as.vector(x)
  }
  if (is.rvobj(x) || is.rvobj(y)) {
    rvmapply(cor, x=x, y=y, ...)
  } else {
    cor(x=x, y=y, ...)
  }
}




