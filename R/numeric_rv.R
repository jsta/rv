

#' @method is.numeric rv
#' @export
is.numeric.rv <- function (x) {
  all(rvsimapply(x, is.numeric))
}

#' @method as.numeric rv
#' @export
as.numeric.rv <- function (x, ...) {
  simapply(x, as.numeric, ...)
}

#' @method as.numeric rvfactor
#' @export
as.numeric.rvfactor <- function (x, ...) {
  simapply(x, as.numeric, ...)
}

