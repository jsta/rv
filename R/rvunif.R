#' Generate Random Vectors from a Uniform Sampling Model
#' 
#' Generates random variables from a Uniform sampling model.
#' 
#' 
#' @param n integer: number of scalars to generate
#' @param min lower limit of the distribution, (may be random)
#' @param max upper limit of the distribution, (may be random)
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples #' 
#'   y <- rvunif(1, min=rvunif(1)-1, rvunif(1)+1) # What marginal distribution does y have now?
#' 
#' @export rvunif
rvunif <- function (n=1, min=0, max=1) {
  rvvapply(runif, n.=n, min=min, max=max)
}

