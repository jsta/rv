.packageName <- "rv"




#' Generate Posterior Simulations for lm or glm Objects (defunct)
#' 
#' 
#' \emph{DEFUNCT.} Use \code{\link{posterior}} instead.
#' 
#' Generate posterior simulations for a given fitted linear or general linear
#' model, assuming the standard "noninformative" priors on the unknowns.
#' 
#' 
#' @aliases postsim postsim.lm postsim.glm
#' @param fit an lm or glm object
#' @return A (named) random vector for each fitted coefficient.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @export postsim
postsim <- function(fit)
{
  warning("'postsim' is defunct. Use 'posterior' instead.")
  UseMethod('postsim')
} 

