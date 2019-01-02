# ========================================================================
# print.rv  -  print summary of a rv on the console
# ========================================================================



#' Print Distribution Summary of a Random Variable
#' 
#' Prints a summary of the random variable object.
#' 
#' Invokes first the summary method of the object, then prints the result.
#' 
#' @param x an rv object
#' @param digits minimal number of significant digits
#' @param \dots further arguments passed to or from other methods
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{rvfactor}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   print(rvnorm(mean=rvnorm(1)))
#' 
print.rv <- function(x, digits=rvpar("print.digits"), ...) {
  if (length(x)==0) {
    return(cat("rv(0)\n"))
  }
  print(summary(x, ...), digits=digits, ...)
}



