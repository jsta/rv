#' Set or Query Parameters of the 'rv' Package
#' 
#' Sets or retrieves parameters of the \code{rv} package.
#' 
#' \describe{ \item{rvcol}{color of a random point (interval), such as 'red' or
#' 'blue'} \item{rvlex}{middle interval expansion factor} \item{rvlwd}{line
#' weight of a random interval} \item{print.digits}{number of digits to show in
#' the summaries} \item{rvpoint}{what to output when plotting a random point;
#' default \code{list("95\%", "50\%", "mean")}} \item{point.sample}{number of
#' points to plot when plotting a rv-rv scatterplot. Default \code{400}.}
#' \item{line.sample}{number of lines to draw when plotting a random sample of
#' lines (see \code{abline}). Default \code{20}.}
#' \item{summary.dimnames}{logical; output dimnames in the summary of an
#' \code{rv} object? Default \code{TRUE}.}
#' \item{summary.quantiles.numeric}{vector of quantiles to compute for the
#' summary of a numeric \code{rv} object.}
#' \item{summary.quantiles.integer}{vector of quantiles to compute for the
#' summary of an integer-valued \code{rv} object. By default contains 0 and 1
#' (for the min and max values).} }
#' 
#' @param \dots arguments in tag = value form, or a list or character vecetor
#' of tagged values. The available tags are described below.
#' @return In the case of a single tag query, the requested value.
#' 
#' In the case of multiple tag query, a list of requested values.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   rvpar()$rvcol
#'   rvpar("rvcol")
#' 
#' @export rvpar
rvpar <- function (...) {
  args <- list(...)
  rvpar  <- getOption("rv")
  Pars <- rvpar$par
  if (length(args)==0) {
    return(Pars)
  }
  arg.names <- names(args)
  if (is.null(arg.names)) {
    args <- unlist(args)
    p <- Pars[args]
    if (length(args)==1) {
      return(p[[args]])
    } else {
      return(p)
    }
  }
  oldpar <- Pars
  for (name in arg.names) { 
    if (nchar(name)>=1) {
      Pars[name] <- args[name]
    }
  }
  rvpar$par <- Pars
  options(rv=rvpar)
  return(oldpar)
}


