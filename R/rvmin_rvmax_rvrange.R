
#' @export
rvmin <- function (x) {
  rvsimapply(x, min, na.rm=FALSE)
}

#' @export
rvmax <- function (x) {
  rvsimapply(x, max, na.rm=FALSE)
}

#' @export
rvrange <- function (x) {
  rvsimapply(x, range, na.rm=TRUE)
}

