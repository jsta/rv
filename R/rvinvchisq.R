#' Generate Random Variables from a Inverse-Chi-Square Sampling Model
#' 
#' \code{rvinvchisq}
#' 
#' \code{rvinvchisq}
#' 
#' @param n integer: number of variables to generate
#' @param df degrees of freedom (may be random)
#' @param scale scale parameter (may be random)
#' @return A random vector (rv object).
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   rvinvchisq(df=3, scale=2)
#' 
#' @export rvinvchisq
rvinvchisq <- function (n=1, df, scale=1) {
  return(scale / (rvchisq(n=n, df=df) / df))
}


