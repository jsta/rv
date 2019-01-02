# ========================================================================
# rvsimapply  -  apply a function to the simulations, componentwise
# ========================================================================



#' Apply a Function to Columns of the Matrix of Simulation of a Random Vector
#' 
#' \code{rvsimapply}
#' 
#' \code{rvsimapply} applies a given function to the \emph{rows} of the
#' simulation matrix of the given random vector.
#' 
#' If the function is to be applied to \emph{rows} of the simulation matrix,
#' use \code{\link{simapply}} or \code{\link{rvmapply}} instead.
#' 
#' Usually used in functions that implement an 'rv'-compatible routine.
#' 
#' @param x an object
#' @param FUN an R function object
#' @param \dots further arguments passed to the function \code{FUN}
#' @return A numeric vector or array.
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
#' @export rvsimapply
rvsimapply <- function(x, FUN, ...)
{
  dx <- dim(x)
  n <- length(x)
  if (n==0) {
    return(NULL)
  }
  mv <- lapply(unclass(x), FUN, ...)
  lmv <- sapply(mv, length)
  if (all(lmv==1)) {
    m <- unlist(mv, use.names=TRUE)
    dim(m) <- dx
    dimnames(m) <- dimnames(x)
    return(m)
  } else if (all(lmv==rvnsims(x))) {
    # simulation-wise function was applied - return an object of same type
    attributes(mv) <- attributes(x)
    return(mv)
  } else if (all(lmv==lmv[1])) {
    m <- unlist(mv)
    m <- matrix(m, nrow=lmv[1], ncol=n)
    dimnames(m) <- list(names(mv[[1]]), names(x))
    return(m)
  } else {
    names(mv) <- names(x)
    return(mv)
  }
}

