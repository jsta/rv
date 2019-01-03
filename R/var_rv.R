# ========================================================================
# var - variance function
# ========================================================================
#

#' Calculate the variance of an rv object
#' 
#' @inheritParams stats::var
#' @param \dots arguments passed to stats::var
#' @export
var.rv <- function(x, ...) {
  simapply(x, stats::var, ...)
}

