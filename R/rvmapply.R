# ========================================================================
# rvmapply  -  apply a function to multiple rv objects
# ========================================================================
# Note. Won't work with functions allowing "blank" arguments
# such as "[" (e.g. x[y,,]). The functions "[" and "[<-" use
# modified versions of simmapply.
#



#' Apply a function to multiple random vector arguments
#' 
#' \code{rvmapply} is the rv-compatible version of \code{\link{mapply}}.  It
#' repeats the function \code{FUN} for each joint draw of the random (or
#' constant) arguments, while allowing vectorizing.
#' 
#' \code{rvmapply} applies a given function to each simulation (vector or
#' array) of the given random vectors, returning a the results as a random
#' vector or array.
#' 
#' The dimensions of each joint draw are preserved.  For an example, see
#' \code{\link{solve}}, that returns the distribution of the inverse of a
#' random matrix.
#' 
#' Usually used in functions that implement an 'rv'-compatible routine.
#' 
#' For an example of a function that uses \code{SAMPLESIZE},
#' \code{\link{abline}}.
#' 
#' @aliases rvmapply rvVectorize
#' @param FUN the function to apply to the simulations of \code{X}.
#' @param MoreArgs Other args passed to \code{FUN} `as is' (must not be rv
#' objects unless the function already accepts them)
#' @param USE.NAMES logical; see \code{\link{mapply}} for details
#' @param SIMPLIFY logical; see \code{\link{mapply}} for details
#' @param SAMPLESIZE if specified, takes a (joint) sample of the simulations
#' and processes only them.
#' @param vectorize.args a character vector of arguments which should be
#' vectorized. Defaults to all arguments to FUN.
#' @param \dots further arguments to \code{FUN}, possibly random vectors or
#' array.
#' @return Depends on \code{FUN}; a random vector or array if \code{FUN} is
#' numeric.
#' @note If the function (\code{FUN}) has an argument ``\code{FUN}", it must be
#' specified within the list supplied to \code{MoreArgs}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{mapply}}, \code{\link{simapply}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords manip
#' @export rvmapply
rvmapply <- function (FUN, ..., MoreArgs=NULL, SIMPLIFY=FALSE, USE.NAMES=TRUE, SAMPLESIZE=NULL) {
  a <- list(...)
  dim.a.names <- dimnames(a)
  a.names <- names(a)
  dimnames(a) <- dim.a.names
  names(a) <- a.names
  a <- lapply(a, .sims.as.list)
  if (!is.null(SAMPLESIZE)) {
    m <- max(sapply(a, length))
    s <- (sample(1:m, size=SAMPLESIZE, replace=TRUE)-1)
    a <- lapply(a, function (x) x[(s %% length(x))+1])
  }
  a <- .Primitive("c")(FUN = FUN, a, SIMPLIFY = FALSE, USE.NAMES=USE.NAMES)
  a$MoreArgs <- MoreArgs
  S <- do.call(mapply, args = a)
  S <- lapply(S, function (x) if (is.null(x)) NA else x) ## DEBUG:: OK??
  r <- rvsims(S)
  ## DEBUG: match the largest-dimensional param in list(...) and 
  ## set the dimnames if they match
  ## if (isTRUE(all.equal(dim(r), dim(x)))) {
  ##  dimnames(r) <- dimnames(x)
  ##}
  return(r)
}



