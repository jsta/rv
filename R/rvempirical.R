#' Generate a Random Vector from an Empirical Distribution
#' 
#' \code{rvempirical} generates a random vector of the same length as data from
#' the empirical distribution of the data.
#' 
#' \code{rvempirical}
#' 
#' @param n Number of i.i.d. rv components to generate
#' @param data Data (constants)
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   y <- c(1.0, 1.2, 3, 1.1, 0.8, 0.9) ## Some data
#'   x <- rvempirical(4, data=y) 
#' 
#' @export rvempirical
rvempirical <- function (n, data) {
  n.sims <- getnsims()
  n.all <- (n.sims * n)
  s <- sample(data, size=n.all, replace=TRUE)
  m <- matrix(s, nrow=n.sims, ncol=n)
  x <- rvsims(m, permute=FALSE)
  return(x)
}
