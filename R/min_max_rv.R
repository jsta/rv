

#' @method min rv
#' @export
min.rv <- function(..., na.rm=FALSE) {
  if (anyisrv(...)) {
    simapply(cbind.rv(...), base::min, na.rm=na.rm)
  } else {
    base::min(..., na.rm=na.rm)
  }
}

#' @method max rv
#' @export
max.rv <- function(..., na.rm=FALSE) {
  if (anyisrv(...)) {
    simapply(cbind.rv(...), base::max, na.rm=na.rm)
  } else {
    base::max(..., na.rm=na.rm)
  }
}

#' @method pmin rv
#' @export
pmin.rv <- function(..., na.rm=FALSE) {
  if (anyisrv(...)) {
    a <- sims(cbind.rv(...), dimensions=TRUE)
    rvsims(t(apply(a, 1, function (m) apply(m, 1, base::pmin))))
  } else {
    base::pmin(..., na.rm=na.rm)
  }
}

#' @method pmax rv
#' @export
pmax.rv <- function(..., na.rm=FALSE) {
  if (anyisrv(...)) {
    a <- sims(cbind.rv(...), dimensions=TRUE)
    rvsims(t(apply(a, 1, function (m) apply(m, 1, base::pmax))))
  } else {
    base::pmax(..., na.rm=na.rm)
  }
}

is.recursive.rv <- function (x) {
  return(FALSE)
}

is.atomic.rv <- function (x) {
  return(TRUE)
}

is.integer.rv <- function (x) {
  return(is.rv(x) && all(rvsimapply(x, is.integer)))
}

is.logical.rv <- function (x) {
  return(is.rv(x) && all(rvsimapply(x, is.logical)))
} 
