#' Generate Random Variables from a Cauchy Sampling Model
#' 
#' Random vector generation for the Cauchy distribution.
#' 
#' For details on the Cauchy distribution, see \link{Cauchy}.  See also
#' \code{\link{rvt}}; Cauchy is a special case of the t-distribution with 1
#' degree of freedom, and therefore \code{rvcauchy(n,location,scale)} is
#' equivalent to \code{rvt(n, mu, scale, df=1)}.
#' 
#' @param n integer: number of variables to generate
#' @param location location parameter (may be random)
#' @param scale scale parameter (may be random)
#' @return A random vector (rv object) of length \code{n}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @export rvcauchy
rvcauchy <- function (n=1, location=0, scale=1) {
  rvvapply(rcauchy, n.=n, location=location, scale=scale)
}


