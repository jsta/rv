#' Componentwise Quantiles of Random Variables
#' 
#' Computes componentwise quantiles of random vectors or arrays.
#' 
#' \code{rvquantile} applies the \code{quantile} function to each column of
#' \code{sims(x)}.
#' 
#' \code{rvmedian} applies \code{median} to the each column of \code{sims(x)}.
#' 
#' @aliases rvquantile rvquantile.rv rvquantile.rvsummary rvmedian
#' @param x an object
#' @param \dots further arguments passed to \code{quantile}
#' @return A \emph{numeric} vector of quantiles.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(3)
#'   rvquantile(x)
#'   rvquantile(x, probs=c(0, 0.01, 0.99, 1))
#'   rvmedian(x)
#' 
#' @export rvquantile
rvquantile <- function(x, ...)
{
  UseMethod("rvquantile")
}

#' @rdname rvquantile
#' @param probs numeric vector of probabilities with values in \emph{[0,1]}
#' @param ignoreInf ignore infinite values
#' 
#' @export
#' @method rvquantile rv
rvquantile.rv <- function(x, probs=c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975), ignoreInf=FALSE, ...)
{
  if (ignoreInf) {
    .f <- function (x) { quantile(x[is.finite(x)], probs=probs, ..., na.rm=TRUE) }
    t(rvsimapply(x, .f))
  } else {
    t(rvsimapply(x, quantile, probs=probs, ..., na.rm=TRUE))
  }
}

#' @method rvquantile rvsummary
#' @export
rvquantile.rvsummary <- function(x, probs=c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975), ...)
{
  Q <- t(sims(x))
  all_probs <- attr(Q, "quantiles")
  M <- NULL
  name <- character(0)
  # if (all(probs %in% all_probs)) ...
  for (p in probs) {
    ix <- (all_probs==p)
    if (any(ix)) {
      M <- cbind(M, Q[,ix,drop=FALSE])
    } else {
      name <- paste(p*100, "%", sep="")
      M <- cbind(M, NA)
      colnames(M)[ncol(M)] <- name
    }
  }
  return(M)
}
