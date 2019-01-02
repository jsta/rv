#' Generate Categorical Random Variables
#' 
#' Generates a random factor (i.e. a categorical random variable), given the
#' probabilities of each category and their corresponding labels.
#' 
#' The length of \code{prob} determines the number of bins.
#' 
#' The vector \code{prob} will be normalized to have sum 1.
#' 
#' @param n integer, number of random variables to generate
#' @param prob vector of probabilities of successes of each trial (may be
#' constant or an rv object)
#' @param levels (character) labels for the categories
#' @return A \emph{random factor} of length \code{length(prob)}.
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
#'   rvcat(1, prob=c(0.5, 0.3, 0.2)) # default levels: 1, 2, 3
#'   rvcat(1, prob=c(5, 3, 2)) # same as above
#'   p <- rvdirichlet(1, alpha=c(0.7, 0.3)) # prior probabilities
#'   rvcat(1, prob=p, levels=c("Group 1", "Group 2"))
#' 
#' @export rvcat
rvcat <- function (n=1, prob, levels=NULL) {
  # NAME
  #  rvcat - Sample Categorical Random Variables
  if (anyisrv(n, prob)) {
    x <- rvmultinom(n=n, size=1, prob=prob)
    s <- sims(x, dimensions=TRUE)
    ds <- dim(s)
    s <- as.logical(s)
    dim(s) <- ds
    f <- function (m) row(m)[m]
    a <- apply(s, 1, f)
    r <- if (is.null(dim(a))) rvsims(a) else rvsims(t(a))
  } else {
    n.sims <- getnsims()
    s <- rmultinom(n=n*n.sims, size=1, prob=prob)
    s <- row(s)[as.logical(s)]
    dim(s) <- c(n.sims, length(s) %/% n.sims)
    r <- rvsims(s)
  }
  rvfactor(r, levels=levels)
}
