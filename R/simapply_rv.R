# ========================================================================
# simapply  -  apply a (numeric) function to the simulations, rowwise, with dimensions
# ========================================================================
# Vectorizes over the simulations of one single rv; 
# for vectorization over a group of rvs, see 'rvmapply', 'rvvapply'.
#



#' Apply a Function to Rows of Simulations of Random Vectors
#' 
#' \code{simapply} applies a given function \code{FUN} to each row of the
#' simulation matrix, returning an rv object.
#' 
#' \code{simapply} applies a given function to the \emph{rows} of the
#' simulation matrix of the given random vector.
#' 
#' If the function accepts \emph{arrays}, use \code{\link{rvmapply}} instead.
#' 
#' If the function is to be applied to each component of the random vector
#' separately (such as in \code{\link{rvmean}}), use \code{\link{rvsimapply}}
#' instead.
#' 
#' Usually used in functions that implement an 'rv'-compatible numeric
#' function.
#' 
#' @param x a random vector.
#' @param FUN a function.
#' @param \dots further arguments passed to \code{FUN}.
#' @return An \code{rv} object, representing the distribution of \code{FUN(x,
#' ...)}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords manip
#' @examples
#' 
#'   #
#'   x <- rvnorm(10)
#'   simapply(x, mean) # Same result as that of mean(x).
#' 
#' @export simapply
simapply <- function(x, FUN, ...) {
  # Works pretty similarly as rvmapply does
  L <- .sims.as.list(x)
  Args <- .Primitive("c")(list(FUN=FUN, SIMPLIFY=FALSE, USE.NAMES=FALSE), list(L))
  Args$MoreArgs <- list(...)
  S <- do.call(mapply, Args)
  r <- rvsims(S) 
  if (isTRUE(all.equal(dim(r), dim(x)))) {
    dimnames(r) <- dimnames(x)
  }
  return(r)
}
