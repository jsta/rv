#' Maxima and Minima of Random Variables
#' 
#' @name Extremes-rv
#' 
#' @description Returns the maxima and minima of the components of a random vector.
#' 
#' \code{rvmin} applies the function \code{min} to each component of the
#' argument \code{x}.  Missing values are removed.
#' 
#' \code{rvmax} applies the function \code{max} to each component of the
#' argument \code{x}.  Missing values are removed.
#' 
#' \code{rvrange} applies the function \code{range} to each component of the
#' argument \code{x}.  Missing values are removed.
#' 
#' \code{min.rv} returns the minimum of the random \emph{vector}, returning
#' thus one random variable. Similarly \code{max.rv} returns the maximum of a
#' vector.
#' 
#' \code{pmin.rv} and \code{pmax.rv} returns the componentwise minima or maxima
#' of several random vectors or constants, yielding thus a random vector of the
#' same length.
#' 
#' @aliases min.rv max.rv pmin.rv pmax.rv
#' @param x an \code{rv} or \code{rvsummary} object
#' @param na.rm remove missing values?
#' @param \dots one or more \code{rv} objects or numeric objects
#' @return A \emph{numeric} vector of the same dimension as \code{x}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{rvmedian}}, \code{\link{rvmean}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvpois(10, lambda=3)
#'   rvmin(x)
#'   rvmax(x)
#'   rvrange(x)
#' 
NULL

#' Get the min values of an rv object
#' @inheritParams Extremes-rv
#' @export
rvmin <- function (x) {
  rvsimapply(x, min, na.rm=FALSE)
}

#' Get the max values of an rv object
#' @export
#' @inheritParams Extremes-rv
rvmax <- function (x) {
  rvsimapply(x, max, na.rm=FALSE)
}

#' Get the value range of an rv object
#' @inheritParams Extremes-rv
#' @export
rvrange <- function (x) {
  rvsimapply(x, range, na.rm=TRUE)
}
