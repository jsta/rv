
#' @export
as.rvobj <- function (x) {
  if (is.rvobj(x)) {
    return(x)
  }
  as.rv(x)
}

#' @export
is.rvobj <- function (x) {
  return(is.rv(x) || is.rvsummary(x))
}

