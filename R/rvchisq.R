#' Generate Random Variables from a Chi-Square Sampling Model
#' 
#' Generates a random vector from a chi-square sampling model.
#' 
#' If any of the arguments are random, the resulting simulations may have
#' non-Poisson marginal distributions.
#' 
#' @param n number of variables to generate
#' @param df integer, degrees of freedom, may be random
#' @param ncp non-centrality parameter, may be random
#' @return A random vector (rv object) of length \code{n}.
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
#' @export rvchisq
rvchisq <- function (n=1, df, ncp = 0) {
  if (missing(ncp)) {
    rvvapply(rchisq, n.=n, df=df)
  } else {
    rvvapply(rchisq, n.=n, df=df, ncp=ncp)
  }
}


