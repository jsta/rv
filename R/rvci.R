#' Credible (Uncertainty) Intervals for Random Scalars
#' 
#' Computes credible (uncertainty) intervals for a given vector, given
#' quantiles or the size of the middle interval
#' 
#' If \code{interval} is of length two or more, the return value will be the
#' quantiles given by \code{range(interval)}.
#' 
#' @param obj random scalar or vector
#' @param interval size of the middle interval or the quantile range of the
#' interval
#' @param one.sided logical, FALSE if two-sided interval is desired
#' @param left logical, indicating if the left one-sided interval is desired
#' @return For two-sided intervals, an array of numbers of dimension
#' \code{c(2,length(x))}, for one-sided intervals, a vector of the same length
#' as \code{x}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   rvci(rvnorm(1), interval=0.683) # Should be about c(-1,1).
#' 
#' @export rvci
rvci <- function(obj, interval=0.95, one.sided=FALSE, left=TRUE)
{
  # NAME
  #  rvci - Uncertainty (Credible) Interval
  # 
  if (length(interval)>1) {
    ci <- rvquantile(obj, probs=range(interval))
  } else if (one.sided) {
    q <- if (left) interval else (1-interval)
    ci <- rvquantile(obj, q)
  } else {
    lower <- (1-interval)/2
    upper <- (lower + interval)
    ci <- rvquantile(obj, c(lower,upper))
  }
  return(ci)
}

