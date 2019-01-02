# ========================================================================
# 
# ========================================================================



#' Random Vectors with a Permutation Distribution
#' 
#' Generates a random vector with each component having a permutation
#' distribution based on the given (fixed) data vector.
#' 
#' 
#' @aliases rvpermut rvpermut
#' @param data a fixed numeric vector
#' @param prob optional probabilities for the components in \code{data}
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvpermut(1:10)
#' 
#' @export rvpermut
rvpermut <- function (data, prob=NULL) {
  ## permutation distribution
  n.sims <- getnsims()
  s <- t(sapply(rep(list(data), n.sims), sample, prob=prob))
  r <- rvsims(s)
  dim(r) <- dim(data)
  names(r) <- names(data)
  return(r)
}

