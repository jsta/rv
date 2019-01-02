
#' @export
rvmedian <- function (x) {
  UseMethod("rvmedian")
}

#' @export
rvmedian.rv <- function (x) {
  rvsimapply(x, median, na.rm=TRUE)
}

#' @export
rvmedian.rvsummary <- function (x) {
  rvquantile(x, probs=0.50)
}
