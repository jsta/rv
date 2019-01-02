#' @name rvmin
#' @noRd
#' @inheritParams Extremes-rv
#' @export
rvmin <- function (x) {
  rvsimapply(x, min, na.rm=FALSE)
}

#' @name rvmax
#' @noRd
#' @export
#' @inheritParams Extremes-rv
rvmax <- function (x) {
  rvsimapply(x, max, na.rm=FALSE)
}

#' @name rvrange
#' @noRd
#' @inheritParams Extremes-rv
#' @export
rvrange <- function (x) {
  rvsimapply(x, range, na.rm=TRUE)
}
