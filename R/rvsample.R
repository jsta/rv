#' Draw a Sample from the Simulation Matrix of a Random Variable
#' 
#' Draws a sample of desired size from each component of a given random
#' variable \code{x}.
#' 
#' Samples (with replacement) from the distribution of the random variable
#' object.  In effect it samples from the rows of the simulation matrix
#' \code{sims(x)}.
#' 
#' @param x an object
#' @param size size of the sample
#' @param jointly return joint simulations and not simulations from each
#' component separately
#' @param reject.na reject each draw that contains an \code{NA}
#' @return A \emph{numeric} array of dimensions \code{size} times
#' \code{length(x)}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   #
#' 
#' @export rvsample
rvsample <- function(x, size=1, jointly=TRUE, reject.na=FALSE)
{
  # NAME
  #   rvsample - Draw Samples from Random Vectors
  #
  xs <- sims(as.rv(x))
  ns <- nrow(xs)
  if (is.null(size) || is.na(size)) size <- ns
  if (jointly) {
    if (reject.na) {
      f <- function (x) any(is.na(x))
      is.na.xs <- apply(xs, 1, f)
      if (all(is.na.xs)) {
        s <- sample(ns, size=size, replace=TRUE, prob=is.na.xs)
      } else {
        s <- sample(ns, size=size, replace=TRUE, prob=!is.na.xs) 
      }
    } else {
      s <- sample(ns, size=size, replace=TRUE)
    }
    s <- xs[s,]
  } else {
    s <- apply(xs, 2, function (s) {
      if (all(nas <- is.na(s))) return(s)
      sample(s[!nas], size=size, replace=TRUE)
    })
  }
  return(s)
}
