	



#' Conditional Random Element Selection
#' 
#' \code{rvifelse} is the \code{rv}-compatible version of the function
#' \code{ifelse}.
#' 
#' \code{rvifelse} returns a \emph{random} value with the same shape as
#' \code{test} which is filled with random or constant elements selected from
#' either \code{yes} or \code{no}, depending on whether the random draw in an
#' element of \code{test} is \code{TRUE} or \code{FALSE}.
#' 
#' @param test an object which can be coerced to logical mode.
#' @param yes return values for true elements of \code{test}
#' @param no return joint simulations and not simulations from each component
#' separately
#' @return A \emph{numeric} array of dimensions \code{size} times
#' \code{length(x)}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{ifelse}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @export rvifelse
rvifelse <- function (test, yes, no) {
  #   rvifelse - If-Then-Else For Random Vectors
  rvmapply(ifelse, test=test, yes=yes, no=no)
}


