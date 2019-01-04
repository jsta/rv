#' Constant Vectors
#' 
#' Functions to coerce or test for non-random objects.
#' 
#' \code{is.constant} returns \code{TRUE} for each component of the argument
#' object if there is only one simulation (that is, the variable is
#' "constant").
#' 
#' Note: rv objects that merely have variance zero are not therefore
#' necessarily "true" constants.
#' 
#' \code{as.constant} coerces \code{rv} or \code{rvsummary} objects into
#' constant strings; \code{NA} is returned for component that is not random.
#' 
#' @return a logical vector (not rv), TRUE if a component is constant w.p. 1
# 
#' @aliases is.constant as.constant as.constant.rv as.constant.rvsummary
#' @param x an object, random variable (rv) or not
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   is.constant(1)         # TRUE
#'   is.constant(as.rv(1))  # TRUE
#'   setnsims(200)
#'   x <- rvbern(prob=0.001)
#'   all(sims(x)==0)        # most probably true
#'   is.constant(x)         # always FALSE
#'   x <- rvnorm(3)
#'   x[1] <- 1
#'   as.constant(x)         # 1, NA, NA
#'   all(is.random(x) & is.na(as.constant(x))) # always TRUE
#' 
#' @export
is.constant <- function(x) {
  # Note: this corresponds to "constant with probability 1", while
  # rvvar(x)==0 would correspond to "constant almost surely"
  return(rvnsims(x)==1)
}

#' Make an rv object have a constant value
#' 
#' @rdname is.constant
#' @export
as.constant <- function(x)
{
  UseMethod('as.constant')
}

#' @method as.constant rv
#' @export
as.constant.rv <- function (x)
{
  z <- rvmean(x)
  z[! is.constant(x)] <- NA
  return(z)
}

as.constant.rvsummary <- function(x)
{
  for (i in which(is.constant(x))) {
    x[[i]] <- x[[i]][1]
  }
  for (i in which(!is.constant(x))) {
    x[[i]] <- NA
  }
  structure(unlist(x), names=names(x))  
}

#' @method as.constant default
#' @export
as.constant.default <- function (x)
{
  return(x)
}
