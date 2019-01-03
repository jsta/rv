## ========================================================================
## function  -  short description
## ========================================================================
## Name:        
## Description: 
## Parameters:  
## Required:    none
## History:     2004-06-  : 
##



#' Generate Posterior Simulations
#' 
#' Generate posterior simulations for a given fitted linear or general linear
#' model, assuming the standard "noninformative" priors on the unknowns.
#' 
#' 
#' @aliases posterior posterior.lm posterior.glm
#' @param obj an object
#' @param \dots further arguments
#' @return A (named) list of random vectors.  For example, the \code{lm} method
#' returns a list with components \code{sigma} (the residual s.d.)  and
#' \code{beta}, the regression coefficients.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   \dontrun{
#'   x <- 1:20
#'   y <- rnorm(length(x), mean=x, sd=10)
#'   print(summary(fit <- lm(y ~ x)))
#'   bayes.estimates <- posterior(fit)
#'   }
#' 
#' @export posterior
posterior <- function(obj, ...) {
  UseMethod("posterior")
} 

## Old code:
##    for (sim in 1:nsim){
##      sigma[sim] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
##      beta[sim,] <- mvrnorm (1, beta.hat, V.beta*sigma[sim]^2)
##    }
##

## ========================================================================
## posterior.lm
## ========================================================================

#' @method posterior lm
#' @export
posterior.lm <- function (obj, ...) {
  ## Modified version of 'sim' (from Andrew Gelman's library (gelman@stat.columbia.edu))
  ##
  summ <- summary(obj)
  sigma.hat <- summ$sigma
  beta.hat <- summ$coefficients[,1]
  V.beta <- summ$cov.unscaled
  k <- summ$df[1]
  n <- k + summ$df[2]
  sigma <- sigma.hat*sqrt((n-k)/rvchisq(1,n-k))
  beta.0 <- as.vector(beta.hat) + sigma * rvnorm(mean=0, var=V.beta)
  names(beta.0) <- names(beta.hat)
  list(beta=beta.0, sigma=sigma)
}

# ========================================================================
# posterior.glm
# ========================================================================

#' @method posterior glm
#' @export
posterior.glm <- function(obj, ...) {
  ## Modified version of 'sim' (from Andrew Gelman's library (gelman@stat.columbia.edu))
  ##
  summ <- summary (obj, correlation=TRUE)
  beta.hat <- summ$coefficients[,1]
  sd.beta <- summ$coefficients[,2]
  corr.beta <- summ$corr
  k <- summ$df[1]
  n <- k + summ$df[2]
  V.beta <- corr.beta * array(sd.beta,c(k,k)) * t(array(sd.beta,c(k,k)))
  ## dimnames(beta) <- list (NULL, names(beta.hat))
  beta <- rvnorm(mean=beta.hat, var=V.beta)
  list(beta=beta, sigma=V.beta)
}

