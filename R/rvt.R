#' Generate Random Variables from a Student-t Sampling Model
#' 
#' Generates a random variable from a Student-t sampling model.
#' 
#' This function generates both univariate (independent and identically
#' distributed) Student-t random variables and multivariate Student-t
#' distributed vectors (with a given scaling matrix).
#' 
#' For details of the parameters, see the entry on \code{mvt} in the
#' \code{mvtnorm} package.
#' 
#' @param n integer, number of scalars to generate
#' @param mu location, may be a rv
#' @param scale scale, may be a rv
#' @param ncp non-centrality parameter
#' @param df degrees of freedom, may be a rv
#' @param Sigma (optional) scaling matrix for multivariate generation
#' @note If any of the arguments are random, the resulting simulations may have
#' non-t marginal distributions.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   df <- 3
#'   x <- rvt(n=1, df=df)
#'   y <- rvnorm(1)/sqrt(rvchisq(1, df=df)/df) # Same distribution as above
#'   print(c(x,y))
#' 
#' @export rvt
#' @importFrom stats rt
rvt <- function (n = 1, mu = 0, scale = 1, df, ncp, Sigma) { ## CHECK
  if (! missing(Sigma)) {
    t <- .rvmvt(n = n, Sigma = Sigma, df = df)
  } else {
    if (missing(ncp)) {
      t <- rvvapply(rt, n. = n, df = df)
    } else {
      t <- rvvapply(rt, n. = n, df = df, ncp=ncp)
    }
    if (scale != 1) 
      t <- (t * scale)
  }
  if (all(mu != 0)) {
    t <- (t + mu) # t + mu, not mu + t (preserves names)
  }
  return(t)
}



.rvmvt <- function (n=1, Sigma, df=1) {
  x <- sqrt(rvchisq(n=n, df=df)/df)
  # DEBUG: But will this work? x is of length n,
  #   but the returned thing is of length n * nrow(Sigma)!
  return(.rvmvnorm(n=n, mean=0, Sigma=Sigma)/x)
}
