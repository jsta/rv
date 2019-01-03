#' Expectation of a Random Variable
#' 
#' \code{rvmean}
#' 
#' \code{rvmean} computes the means of the simulations of all individual
#' components of a random vector (rv) object.
#' 
#' \code{E} is an alias for \code{rvmean}, standing for ``Expectation.''
#' 
#' \code{Pr} is another alias for \code{rvmean}, standing for ``Probability
#' of''; suggested to be used when the argument is a logical statement
#' involving random variables (that is, a description of an event such as
#' \code{x>0} or \code{x>y}). Then \code{Pr(x>0)} gives the probability of the
#' event ``x>0''. The statement \code{x>0} returns a Bernoulli (indicator)
#' random variable object (having 1/0 or TRUE/FALSE values) and the expectation
#' of such variable is just the probability of the event where the indicator is
#' one.
#' 
#' @aliases rvmean E Pr
#' @param x an rv object
#' @return A \emph{numerical} vector with the same dimension as \code{x}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{mean.rv}}: distribution of the arithmetic mean of a
#' vector; \code{\link{rvmin}}, \code{\link{rvmax}}, \code{\link{rvmedian}},
#' \code{link{rvvar}}, \code{\link{rvsd}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(mean=(1:10)/5, sd=1)
#'   rvmean(x)  # means of the 10 components
#'   E(x)       # same as rvmean(x)
#'   Pr(x>1)    # probabilities that each component is >1.
#' 
#' @export
rvmean <- function (x) {
  UseMethod("rvmean", x)
}

#' @method rvmean rv
#' @export
rvmean.rv <- function (x) {
  m <- colMeans(sims(x), na.rm=TRUE)
  names(m) <- names(x)
  dim(m) <- dim(x)
  dimnames(m) <- dimnames(x)
  return(m)
}

#' @method rvmean rvsummary
#' @export
rvmean.rvsummary <- function (x) {
  unlist(rvattr(x, "mean"), use.names=TRUE)
}

#' @method rvmean default
#' @export
rvmean.default <- function (x) {
  if (!is.numeric(x)) {
    x[] <- as.numeric(x)
  }
  return(x)
}

#' @export
E <- function (x) {
  rvmean(x)
}

#' @noRd
#' @export
#' @param X a logical rv object
Pr <- function (X) {
  if (! is.logical.rv(X)) {
      stop("Argument for Pr must be a logical statement such as 'x>0'")
  }
  rvmean(X)
}
