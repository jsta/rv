#' Detach the rv package
#' 
#' \code{detachrv} detaches the rv package and restores the original functions
#' in \code{base}, \code{graphics} and \code{stats} packages.
#' 
#' Currently equivalent to \code{detach("package:rv")}.
#' 
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   \dontrun{
#'    library(rv)
#'    detachrv()
#'   }
#' 
#' @export detachrv
detachrv <- function ()
{
  detach("package:rv")
  unloadNamespace("rv")
}

