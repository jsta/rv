
#' @export
rvmedian <- function (x) {
  UseMethod("rvmedian")
}

#' @export
#' @method rvmedian rv
rvmedian.rv <- function (x) {
  rvsimapply(x, median, na.rm=TRUE)
}

#' @export
#' @method rvmedian rvsummary
rvmedian.rvsummary <- function (x) {
  rvquantile(x, probs=0.50)
}
