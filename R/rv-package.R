#' Distributions of Various Statistics of a Random Vector and Array
#' 
#' \code{var.rv}, \code{cov.rv} and \code{cor.rv} compute the distribution of
#' the variance statistic of x and the distribution of the covariance statistic
#' or the correlation statistic of x and y if these are vectors.  If x and y
#' are matrices then the covariances (or correlations) between the columns of x
#' and the columns of y are computed.
#' 
#' 
#' These functions are compatible with \emph{both} numeric and rv objects.  To
#' make your code compatible with \code{rv} objects, use e.g. \code{sd.rv}
#' instead of \code{sd}.
#' 
#' The functions \code{cov.rv} is implemented by applying the corresponding
#' numerical function to the rows of the simulation matrices of \code{x} and
#' \code{y} and forming a new \code{rv} object from the resulting vector of
#' simulations.  Alternatively \code{x} may be a random matrix (and \code{y}
#' \code{NULL}).  %Then the numerical function \code{cov}.
#' 
#' \code{cor.rv} works similarly, but returns the distribution of the
#' correlation statistic (i.e. function).
#' 
#' \code{var.rv} computes the distribution of the variance statistic.
#' \code{sd.rv} is the square root of the result obtained by \code{var.rv}.
#' 
#' @aliases cov.rv cor.rv var.rv sd.rv
#' @param x a numeric or random vector, matrix, or a data frame
#' @param y \code{NULL} (default) or a vector, matrix or data frame with
#' compatible dimensions to x. The default is equivalent to y = x (but more
#' efficient).
#' @param \dots further arguments passed to the corresponding numeric functions
#' @name distrib_rv
#' @return A random vector or array.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords internal
#' @examples
#' 
#'   #
#' 
NULL

#' Simulation-based Random Variable Objects
#' 
#' `\code{rv}' implements a simulation-based random variable object class.
#' 
#' Please refer to the vignette: \code{vignette("rv")} for details.
#' 
#' \tabular{ll}{ Package: \tab rv \cr Version: \tab 2.3.0 \cr Date: \tab
#' 2013-05-18 \cr Namespace: \tab rv \cr Depends: \tab R(>= 2.10.0), methods,
#' utils, grDevices, graphics \cr License: \tab GPL-2 \cr }
#' 
#' @name rv-package
#' @docType package
#' @author Jouni Kerman \email{jouni@@kerman.com} Package built on Sat May 18
#' 22:47:25 CEST 2013
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
NULL
