# ========================================================================
# rvsims  -  generate a vector of rvs from a simulation matrix
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
# If the sims array is 1-dimensional,
# it is taken to be the vector of simulations for one variable;
# if sims is 2-dimensional (L x n)
# it contains L simulations for each of n variables.
# if the sims argument is a list, assumes that each component of the
# list is a single draw from a distribution. A component may be a list,
# when the list will be recursively built into a list of rv objects.
#

.recycle_index <- function (length.from, length.to) {
  if (length.from>=1 && length.to>=1) {
    return(1+(0:(length.to-1)) %% length.from)
  }
  return(0)
}

resizeSims <- function (s, vector.length, n.sims) # NOEXPORT
{
  # here the simulations must be in columns
  if (!missing(vector.length)) {
    nr <- nrow(s) # params
    if (nr!=vector.length) {
      ix <- .recycle_index(nr, vector.length)
      s <- s[ix, , drop=FALSE]
    }
  }
  nc <- ncol(s) # sims
  if (nc!=1 && !missing(n.sims)) {
    if (nc!=n.sims) {
      ix <- .recycle_index(nc, n.sims)
      s <- s[, ix, drop=FALSE]
    }
  }
  return(s)
}

.rvsims.list <- function (x, n.sims=getnsims(), permute = FALSE) 
{
  # Assume that all elements in the list have the same dimensions.
  # This may be modified later -- filling with NAs
  dx <- dim(x[[1]])
  s <- sapply(x, length)
  if (!all(s==s[1])) {
    # some elements had different dimensions!
    stop("Simulation list was not consistent")
  }
  S <- t(matrix(unlist(x, recursive=FALSE), nrow=s[1]))
  if (!is.list(S)) {
    r <- rvsims(S, n.sims=n.sims, permute=permute)
    dim(r) <- dx
    names(r) <- names(x[[1]])
    return(r)
  }
  r <- list()
  for (i in 1:ncol(S)) {
    r[[i]] <- .rvsims.list(S[,i])
  }
  names(r) <- names(x[[1]])
  return(r)
}



#' Create Random Vectors from Simulation Draws
#' 
#' \code{rvsims} takes a vector, matrix, or list (\code{sims}) containing
#' simulations, and returns a random vector (an object of type 'rv')
#' 
#' If \code{sims} is a plain numeric vector, this is interpreted to be
#' equivalent to a one-dimensional array, containing simulations for one single
#' random variable.
#' 
#' If the array \code{sims} is one-dimensional, this is interpreted to be
#' equivalent to a two-dimensional array with 1 column.
#' 
#' If \code{sims} is two-dimensional, the \emph{columns} are supposed to
#' contain simulations for one or more several random variables.
#' 
#' If \code{sims} is a list, the numeric vectors are recursively combined to a
#' list of random vectors: each component of the list is supposed to be
#' containing \emph{one} (joint) draw from some distribution---this may be a
#' list.
#' 
#' If \code{permute} is \code{TRUE}, the simulations are scrambled, i.e. the
#' joint draws are permuted randomly.
#' 
#' @param sims an array of simulations (1, or 2-dimensional) or a list
#' @param n.sims number of simulations to save
#' @param permute logical, indicate if scramble the simulations
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   ## x and y have the same distributions but not the same simulations:
#'   n.sims <- 200L
#'   setnsims(n.sims)
#'   y <- rvnorm(1)
#'   x1 <- rvsims(rnorm(n.sims))
#'   ##
#'   s <- sims(x1)
#'   z <- array(s) ## One-dimensional array
#'   x2 <- rvsims(z) ## Same as 
#'   ##
#'   identical(x1, x2) ## TRUE
#'   ##
#'   s <- t(array(rnorm(n.sims * 2, mean=c(0, 10)), dim=c(2, n.sims)))
#'   x3 <- rvsims(s)
#'   identical(2L, length(x3)) ## TRUE
#' 
#' @export rvsims
rvsims <- function(sims, n.sims = getnsims(), permute = FALSE) {
  if (is.list(sims)) {
    if (is.object(sims)) {
      stop(sprintf("rvsims: Cannot apply 'rvsims' on an object of class '%s'", class(sims)[1]))
    }
    return(.rvsims.list(sims, n.sims = n.sims, permute = permute))
  }
  is_factor <- (is.character(sims) || is.factor(sims))
  if (is.factor(sims)) {
    sims[TRUE] <- as.character(sims)
  }
  if (length(sims) < 1) {
    return(rv(0))
  }
  d.s <- dim(sims)
  if (length(d.s) >= 3L) {
    n.sims    <- d.s[1]
    dim.rest  <- d.s[-1]
    n.rest    <- prod(dim.rest)
    dim(sims) <- c(n.sims, n.rest)
    d.s       <- dim(sims)
  } else {
    dim.rest <- integer(0)
  }
  if (! length(d.s) %in% 0:2) {
    stop("rvsims: Argument must be a vector or an array or a list")
  }
  if (is.null(d.s <- dim(sims)) || length(d.s) == 1L) {
    d.s <- dim(sims) <- c(length(sims), 1)
  }
  n.sims.max <- d.s[1]
  if (length(d.s)==2) { # A regular matrix of simulations
    .names <- dimnames(sims)[[2]]
  }
  if (n.sims.max > 1) {
    if  (nrow(sims)!=n.sims) {
      sims <- t(resizeSims(t(sims), n.sims=n.sims))
    }
    if (permute) {
      .order <- sample(n.sims)
      sims <- sims[.order, , drop=FALSE]
    }
  }
  vec <- split(sims, col(sims))
  names(vec) <- .names
  class(vec) <- class(rv())
  if (length(dim.rest) > 0) {
    dim(vec) <- dim.rest
  }
  if (is_factor) {
    as.rvfactor(vec)
  } else {
    return(vec)
  }
}
