#' Numeric Random Vectors
#' 
#' Creates or coerces rv objects of type "numeric".
#' 
#' \code{is.numeric(x)} returns \code{TRUE} if and only if \emph{each}
#' component of \code{x} is numeric-valued.
#' 
#' \code{as.numeric.rv} coerces an rv object into numeric-valued one.  In
#' effect, the function \code{as.numeric} is applied to all simulations.
#' 
#' Random factors are not numeric (just as non-random factors aren't).
#' 
#' @aliases numeric.rv as.numeric.rv is.numeric.rv is.numeric.rvfactor
#' as.numeric.rvfactor
#' @param x an rv object to be coerced or tested.
#' @param \dots further arguments passed to or from other methods.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @name numeric_rv
#' @seealso \code{\link{numeric}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- as.logical(rvbern(1,0.5)) # Bernoulli rv
#'   is.numeric(x)           # FALSE
#'   x <- as.numeric(x)      # coerce to numeric; all TRUEs become ones, FALSEs zeros
#'   is.numeric(x)           # TRUE
#' 
NULL

#' @method is.numeric rv
#' @rdname numeric_rv
#' @export
is.numeric.rv <- function (x) {
  all(rvsimapply(x, is.numeric))
}

#' @method as.numeric rv
#' @rdname numeric_rv
#' @export
as.numeric.rv <- function (x, ...) {
  simapply(x, as.numeric, ...)
}

#' @method as.numeric rvfactor
#' @rdname numeric_rv
#' @export
as.numeric.rvfactor <- function (x, ...) {
  simapply(x, as.numeric, ...)
}
