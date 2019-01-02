#' Distribution of Order Statistics of a Random Vector
#' 
#' \code{sort.rv} computes the distribution of the order statistics of a random
#' vector.
#' 
#' The result is the \emph{distribution} of the order statistic of the given
#' vector \code{x}: that is, the \code{sort} function is applied to each
#' \emph{row} of the matrix of simulations of \code{x} (\code{sims(x)}) and
#' returned then in random vector form.
#' 
#' See \code{\link{sort}} for further details of the function \code{sort}.
#' 
#' @param x a random vector
#' @param \dots further arguments passed to \code{sort.rv}
#' @return An rv object of the same length as \code{x}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{sort}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords manip
#' @examples
#' 
#'   #
#' 
sort.rv <- function (x, ...) {
  simapply(x, base::sort, ...)
}

