#' Variances of Components of Random Vectors
#' 
#' Computes variances of the simulations of components of a random vector of
#' array.
#' 
#' \code{rvvar} computes the means of the simulations of all individual
#' components of a random vector (rv) object.
#' 
#' That is, \code{rvvar} applies the function \code{var} to the vector of
#' simulations of each component of \code{x}, thus computing "columnwise"
#' variances of the matrix of simulations of \code{x}.
#' 
#' \code{rvsd} applies the function \code{sd} to the vector of simulations of
#' each component of \code{x}, thus computing "columnwise" standard deviations
#' of the matrix of simulations of \code{x}.
#' 
#' @aliases rvvar rvsd rvsd.rv rvsd.rvsummary rvvar.rv rvvar.rvsummary
#' @param x an object
#' @return A numeric vector or array (of the same dimension as that of
#' \code{x})
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{rvmin}}, \code{\link{rvmax}}, \code{\link{rvmedian}},
#' \code{\link{rvsd}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(mean=0, var=1:10)
#'   rvvar(x)
#'   rvsd(x)
#' 
#' @export rvvar
rvvar <- function (x)
{
  UseMethod("rvvar")
}

rvvar.rv <- function (x) # NOEXPORT
{
  S <- sims(x)
  m <- colMeans(S, na.rm=TRUE)
  ns <- rvnsims(x)
  v <- ((colSums(S^2)-ns*(m^2))/(ns-1))
  v[ns==1] <- 0
  names(v) <- names(x)
  dim(v) <- dim(x)
  dimnames(v) <- dimnames(x)
  return(v)
}

rvvar.rvsummary <- function (x) # NOEXPORT
{
  return(unlist(rvattr(x, "sd"), use.names=TRUE)^2)
}

rvvar.default <- function (x) # NOEXPORT
{
  rep.int(0, length(x))
}

rvsd <- function (x)
{
  UseMethod("rvsd")
}

rvsd.rv <- function (x) # NOEXPORT
{
  sqrt(rvvar(x))
}

rvsd.rvsummary <- function (x) # NOEXPORT
{
  unlist(rvattr(x, "sd"), use.names=TRUE)
}

rvsd.default <- function (x) # NOEXPORT
{
  rep.int(0, length(x))
}



