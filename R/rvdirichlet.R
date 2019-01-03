#' Generate Random Variables from a Dirichlet Sampling Model
#' 
#' Generates random variables from a Dirichlet sampling model.
#' 
#' The Dirichlet distribution is a generalization of the Beta distribution.
#' (If alpha is of length two, \code{rvdirichlet} draws from the Beta model.)
#' 
#' @param n integer: number of vectors to generate
#' @param alpha the parameter vector; may be random
#' @return A random vector (rv object) of length \code{n}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples \dontrun{
#' 
#'   a <- rvdirichlet(1, alpha=c(6, 3, 1)) # 
#'   sum(a) # one with probability 1   
#'   }
#' 
#' @export rvdirichlet
rvdirichlet <- function (n = 1, alpha)  {
  x <- NULL
  for (i in 1:n) {
    g <- rvgamma(n = 1, shape = alpha, scale = 1)
    x <- cbind.rv(x, g/sum(g))
  }
  return(x)
}
