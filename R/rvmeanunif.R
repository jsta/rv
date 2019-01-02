#' The distribution of the mean of uniform random variables
#' 
#' The distribution of the mean of uniform random variables with each of them
#' in the interval \code{(-1, 1)}, then scaled and shifted.
#' 
#' Assuming that all inputs are constants, each generated variable has a mode
#' (center) at \code{mode}, constrained between \code{(-scale, scale)}.
#' 
#' The shape becomes more and more bell-shaped (Normal) as the number of the
#' independent variables in the sum (mean) increases.
#' 
#' The case of \code{df=2} (mean of two variables) is the special case of the
#' symmetric triangular distribution in the range
#' 
#' @aliases rvmeanunif rvtriang
#' @param n Length of the vector to output
#' @param mode Mode (center) of the distribution
#' @param scale Scale (half-width) of the distribution around the mode
#' @param df ``degrees of freedom'': number of independent components to
#' average
#' @return A random vector of length \code{n}.
#' @author J Kerman
#' @keywords dist
#' @examples
#' 
#'   x <- rvtriang(1)
#'   y <- rvmeanunif(df=2) ## same distribution as that of x
#' 
#' @export rvmeanunif
#' @importFrom stats runif
rvmeanunif <- function (n=1, mode=0, scale=1, df) {
  if (n == 1) {
    u <- rvunif(n=df, min=-1, max=1)
    x <- mean(u)
  } else if (n > 1) {
    n.sims <- getnsims()
    M <- array(runif(n=n.sims * n * df, min=-1, max=1), c(n.sims, n, df))
    R <- t(apply(M, MARGIN=1, rowMeans))
    x <- rvsims(R)
  } else {
    stop("n<1")
  }
  return(mode + scale * x)
}

rvtriang <- function (n=1, mode=0, scale=1) {
  rvmeanunif(n=n, mode=mode, scale=scale, df=2)
}

